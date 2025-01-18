import openai
import json
import os
import random
import pandas as pd
import re
import requests
from datetime import datetime
import requests
from xml.etree import ElementTree as ET

from utils import mkdirsafe
from utils import get_gene_summary_citations, get_family_summary_citations
from utils import find_citations

from utils import get_selected_genes_filepaths

# needed for specificity filtering
# from utils_famfilter import get_cross_references, get_stats_for_gene_name

model_cost = {'gpt-4o-2024-05-13': 5,
              'gpt-4-turbo-2024-04-09': 10
             }


def init_gpu_usage():
    return dict({'NUM_CALLS': 0, 'NUM_CHARS': 0})


def verbose_gpu_usage(GPT_USAGE, config):
    cost_per_1M = model_cost[config['model']]
    print('GPT-4 usage:')
    # 4 chars ≈ 1 token
    print(f'{GPT_USAGE["NUM_CALLS"]} calls, {GPT_USAGE["NUM_CHARS"]} tokens processed (tokens ≈ characters / 4). Approximate cost: ${round(GPT_USAGE["NUM_CHARS"] / 4 / 1000000 * cost_per_1M, 2)}')


def create_summary(client, messages, model):
    try:
        completion = client.chat.completions.create(
            messages=messages,
            model=model,
        )
    except openai.NotFoundError:
        raise
    except openai.BadRequestError:
        raise
    except openai.OpenAIError:
        raise
    else:
        return completion.model_dump(exclude_unset=True)

def parse_response(response):
    choice = response["choices"][0]
    if choice["finish_reason"] == "stop":
        message = choice["message"]["content"]
        return message
    else:
        return 'No response'


def get_gpt_genes_response(client, query, dir_name, gpt4_n, config, prompt_th):
    print('INFO: running per-gene summarization')
    model = config['model']
    GPT_USAGE = init_gpu_usage()
    
    prefix = f'{dir_name}/{query}/per_gene_joined_prompts_dirty'
    response_df_path = f'{dir_name}/{query}/per_gene_summaries_{model}_sample_{gpt4_n}_4o.csv'

    # maybe i can make another version of the response instead
    if os.path.exists(response_df_path):
        print(f'INFO: Gene summaries exist for this sample size ({gpt4_n}) and model ({model})! No rerunning GPT-4 will happen.')
        df = pd.read_csv(response_df_path)
        get_gene_summary_citations(query, dir_name, response_df_path)
        return GPT_USAGE, list(df['gene_name']), list(df['response'])

    # SHOULD I MOVE IT?
    gene_prompt_paths = get_selected_genes_filepaths(query, dir_name, gpt4_n)
    gene_names, gpt4_responses = [], []
    
    if len(gene_prompt_paths) == 0:
        print('No genes found!')
    else:
        for query_path in gene_prompt_paths:
            df = pd.read_csv(prefix + '/' + query_path)
            prompt = list(df['prompt'])[0]
            num_snippets = int(list(df['num_snippets'])[0])
            print(query_path, 'length of prompt:', len(prompt), 'chars, num snippets:', num_snippets)
    
            messages=[
                {"role": "system", "content": config['system']},
                {"role": "user", "content": prompt}
            ]

            print('INFO: Running GPT4 API for gene summarization')
    
            response = create_summary(client, messages, model)
            GPT_USAGE['NUM_CALLS'] += 1
            GPT_USAGE['NUM_CHARS'] += len(prompt)
            parsed_response = parse_response(response)

            gene_name = os.path.splitext(query_path)[0]
            gene_names.append(gene_name)
            
            gpt4_responses.append(parsed_response)
    
        df = pd.DataFrame({'gene_name': gene_names,
                           'response': gpt4_responses})
        
        df.to_csv(response_df_path)
        
    get_gene_summary_citations(query, dir_name, response_df_path)
    
    return GPT_USAGE, gene_names, gpt4_responses
    

def make_family_prompt(gene_names, gpt4_responses, config):
    # prompt_text = config['instr_1_fam_temp'].format(len(gene_names)) + '\n' + '[CONT]' + '\n'
    prompt_text = config['instr_1_fam_temp'].format(len(gene_names)) + '\n' + '\n'
    
    for i in range(len(gene_names)):
        gene_name, gene_response = gene_names[i], gpt4_responses[i]
        prompt_text += f'Gene {i + 1}: {gene_name}\n'
        prompt_text += gene_response 
        prompt_text += '\n'
    
    # prompt_text += '[/CONT]\n'
    prompt_text += '\n'
    if config['instr_2_fam_temp'].find('file:') != -1:
        with open(config['instr_2_fam_temp'][5:], 'r') as f:
            instr_2_fam_temp = f.read()
    else:
        instr_2_fam_temp = config['instr_2_fam_temp']
    prompt_text += instr_2_fam_temp
    
    return prompt_text


BAD_SNIPPET_ERROR_STRING = 'PMID:bad_snippet'
NCBI_PARSE_ERROR_STRING = 'PMID:error_parsing_NCBI_response'
NCBI_FETCH_ERROR_STRING = "PMID:error_fetching_data_from_NCBI"
ERROR_STRINGS = ('[' + BAD_SNIPPET_ERROR_STRING + ']',
                 '[' + NCBI_PARSE_ERROR_STRING + ']',
                 '[' + NCBI_FETCH_ERROR_STRING + ']')


def pmcid_to_pmid(pmc_id):
    service_root = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
    email = "iponamareva@ebi.ac.uk"
    url = service_root + f'?email={email}&ids={pmc_id}'
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            root = ET.fromstring(response.content)
            pmid = root.find("record").attrib['pmid']
            pmid = 'PMID:' + pmid
        except:
            pmid = NCBI_PARSE_ERROR_STRING
    else:
        pmid = NCBI_FETCH_ERROR_STRING
    return pmid


def normalize_matches(matches, error_strings):
    d = dict()
    for match in matches:
        s = set(match.split(', '))
        for error_string in error_strings:
            if error_string in s:
                s.remove(error_string)
                
        temp = ', '.join(list(s))
        d[match] = temp
    return d
        

def substitute_pmcid_to_pmid(text, df):
    if len(df) == 0:
        print('WARNING: No citations found. This can be an error.')
        return text
        
    citations = set(find_citations(text))
    snippet_error_status =  dict(zip(df['snippet_id'], df['snippet_error']))
    
    for snippet_id in citations:
        idx = text.find(snippet_id)
        if snippet_error_status[snippet_id] == False:
            pmcid = snippet_id.split('_')[0]
            # can save them
            pmid = pmcid_to_pmid(pmcid)
            text = text.replace(snippet_id, pmid)
        else:
            text = text.replace(snippet_id, BAD_SNIPPET_ERROR_STRING)

    pattern = r"\[PMID:\w+\](?:, \[PMID:\w+\])*"
    matches = re.findall(pattern, text)

    d_normed_matches = normalize_matches(matches, error_strings=ERROR_STRINGS)
    for key in d_normed_matches:
        text = text.replace(key, d_normed_matches[key])
        
    return text


def save_response(response, dir_name, query, model, N, timestamp_str, type='raw', response_dir='ALL_RESPONSES'):
    response_txt_path = f'{dir_name}/{query}/RESPONSE_{type}_{model}_sample_{N}_{timestamp_str}.txt'
    with open(response_txt_path, 'w') as f:
        print(response, file=f)

    mkdirsafe(response_dir)
    
    response_txt_path = f'{response_dir}/{dir_name}_{query}_RESPONSE_{type}_{model}_sample_{N}_{timestamp_str}.txt'
    with open(response_txt_path, 'w') as f:
        print(response, file=f)


def get_gpt_family_response(client, query, dir_name, gene_names, gpt4_responses, config):
    print('INFO: Running family summarization')
    # allowing multiple responses
    
    model = config['model']
    GPT_USAGE = init_gpu_usage()
    
    N = len(gene_names)
    
    now = datetime.now()
    timestamp_str = now.strftime("%Y-%m-%d %H:%M:%S")
        
    prompt_text = make_family_prompt(gene_names, gpt4_responses, config)
    with open(f'{dir_name}/{query}/PROMPT_family_text_sample_{N}_{timestamp_str}.txt', 'w') as f:
        print(prompt_text, file=f)

    parsed_response, pmid_parsed_response = 'No response', 'No response'
        
    messages=[
        {"role": "system", "content": config['system']},
        {"role": "user", "content": prompt_text}
    ]

    response = create_summary(client, messages, model)
    GPT_USAGE['NUM_CALLS'] += 1
    GPT_USAGE['NUM_CHARS'] += len(prompt_text)
    parsed_response = parse_response(response)

    save_response(parsed_response, dir_name, query, model, N, timestamp_str, type='raw')

    df_factcheck = get_family_summary_citations(query, dir_name, parsed_response)
    pmid_parsed_response = substitute_pmcid_to_pmid(parsed_response, df_factcheck)
    
    save_response(pmid_parsed_response, dir_name, query, model, N, timestamp_str, type='pmid')
    
    return GPT_USAGE, parsed_response, pmid_parsed_response


# not used now
def factcheck_gene_summary(client, query, dir_name, gpt4_n, config):
    snippet_ids, snippet_texts, factcheck_results, summaries_for_df = [], [], [], []
    
    model = config['model']
    GPT_USAGE = init_gpu_usage()
    response_df_path = f'{dir_name}/{query}/per_gene_summaries_{model}_sample_{gpt4_n}_4o.csv'
    df_summaries = pd.read_csv(response_df_path)
    data = get_gene_summary_citations(query, dir_name, response_df_path)

    gene_names, summaries = list(df_summaries['gene_name']), list(df_summaries['response'])

    for i in range(len(gene_names)):
        gene_name, summary = gene_names[i], summaries[i]
        print(gene_name)
        for snippet_id in data[gene_name]:
            prompt = config['factcheck'] + '\n'
            prompt += f'Text 1: {summary}\n'
            prompt += f'Text 2: {data[gene_name][snippet_id]}\n'

            messages=[
                    {"role": "system", "content": config['system']},
                    {"role": "user", "content": prompt}
                ]
    
            print('Running GPT4 API')
    
            # response = create_summary(client, messages, model)
            GPT_USAGE['NUM_CALLS'] += 1
            GPT_USAGE['NUM_CHARS'] += len(prompt)
            # parsed_response = parse_response(response)
            parsed_response = 'lorem ipsum'
            
            snippet_ids.append(snippet_id)
            snippet_texts.append(data[gene_name][snippet_id])
            factcheck_results.append(parsed_response)
            summaries_for_df.append(summary)
    
    df = pd.DataFrame({'summary': summaries_for_df,
                       'snippet_id': snippet_ids,
                       'snippet': snippet_texts,
                       'factcheck_result': factcheck_results
                      })
    df.to_csv(f'{dir_name}/{query}/factcheck_gene_results.csv')

    return GPT_USAGE
            

def factcheck_summary(client, query, dir_name, text, config):
    snippet_ids, snippet_texts, factcheck_results = [], [], []
    model = config['model']
    GPT_USAGE = init_gpu_usage()
    
    data = get_family_summary_citations(query, dir_name, text)

    for snippet_id in data:
        prompt = config['factcheck'] + '\n'
        prompt += f'Text 1: {text}\n'
        prompt += f'Text 2: {data[snippet_id]}\n'

        messages=[
                {"role": "system", "content": config['system']},
                {"role": "user", "content": prompt}
            ]

        # print(messages)
        print('Running GPT4 API')

        response = create_summary(client, messages, model)
        GPT_USAGE['NUM_CALLS'] += 1
        GPT_USAGE['NUM_CHARS'] += len(prompt)
        parsed_response = parse_response(response)
        snippet_ids.append(snippet_id)
        snippet_texts.append(data[snippet_id])
        factcheck_results.append(parsed_response)

    df = pd.DataFrame({'summary': [text] * len(snippet_ids),
                       'snippet_id': snippet_ids,
                       'snippet': snippet_texts,
                       'factcheck_result': factcheck_results
                      })
    df.to_csv(f'{dir_name}/{query}/factcheck_results.csv')

    return GPT_USAGE
        