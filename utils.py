import pandas as pd
import re
import requests
from requests.adapters import HTTPAdapter, Retry

import random
import json
import os
import argparse
import pickle
import xml.etree.ElementTree as ET

from collections import defaultdict

from utils_famfilter import make_spec_stats_file

retries = Retry(total=5,
                backoff_factor=0.1,
                status_forcelist=[408, 500, 501, 502, 503, 504])


def get_dbname_from_acc(accession):
    if accession[:4] == 'PTHR': return 'panther'
    if accession[:2] == 'PF': return 'pfam'
    if accession[:3] == 'SSF': return 'ssf'
    if accession[:2] == 'cd': return 'cdd'
    if accession[:5] == 'G3DSA': return 'cathgene3d'


def make_and_save_data(name, output_path, found_mentions):    
    paper_nums, paper_ids, snippet_nums = [], [], []
    snippets = []
    for i, paper_id in enumerate(found_mentions):
        for j, entry in enumerate(found_mentions[paper_id]):
            paper_nums.append(i)
            paper_ids.append(paper_id)
            snippet_nums.append(j)
            snippets.append(entry)
            
    df_domain_mentions = pd.DataFrame({
        'paper_id': paper_ids,
        'paper_num': paper_nums,
        'snippet_num': snippet_nums,
        'snippet': snippets,
    })
    
    df_domain_mentions.to_csv(output_path)
    print(f'LOG: Saved snippets: {output_path}, found {len(snippets)} snippets')


def mkdirsafe(path):
    if not os.path.exists(path):
        os.mkdir(path)


def format_name(domain_name):
    tmp = '_'.join(domain_name.split())
    tmp = tmp.replace('/', '_')
    return tmp


def verbose_args(args):
    print('-'*80)
    for arg in vars(args):
        print(arg, getattr(args, arg))
    print('-'*80)


def save_args_log(args, config, run_name):
    with open(f'{args.dir_name}/{args.query}/{run_name}/run_args.log', 'w') as f:
        for arg in vars(args):
            print(arg, getattr(args, arg), file=f)
        print('CONFIG', file=f)
        for key in config:
            print(key, config[key], file=f)


def get_first_page(query_text, sess):
    r = sess.get(f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="{query_text}"&format=json&resultType=core')
    
    if r.status_code != 200:
        print('ERROR:', r.status_code, flush=True)
        print('Currently no handling happens')
    
    # error can happen here. todo: handler.
    try:
        data = r.json()
    except requests.exceptions.JSONDecodeError:
        return None
    
    return data


def get_intervals(idxs, window):
    if len(idxs) == 0:
        return []
        
    coords = []
    cur_start_coord = -1
    cur_end_coord = -1
    for idx in idxs:
        idx_start = idx - window
        idx_end = idx + window
        
        # if we find overlap, extend interval
        if idx_start <= cur_end_coord:
            cur_end_coord = idx_end
        # else, save previous interval and create new one
        else:
            if (cur_start_coord != -1 and cur_end_coord != -1):
                coords.append((cur_start_coord, cur_end_coord))
            cur_start_coord = idx_start
            cur_end_coord = idx_end
            
    if (cur_start_coord != -1 and cur_end_coord != -1):
        coords.append((cur_start_coord, cur_end_coord))
    return(coords)


def parse_xml_response(domain_name, text, window):
    try:
        root = ET.fromstring(text)
    except:
        print('INFO: XML tree parsing error in gene name', domain_name)
        return []
        # raise xml.etree.ElementTree.ParseError

    all_text = ''
    selected_text = []
    
    for child in root.iter('p'):
        text_content = ''.join(child.itertext())
        all_text += text_content
        all_text += '\n'

    try:
        idxs = [m.start() for m in re.finditer(domain_name, all_text, flags=re.I)]
    except:
        print('INFO: Regex error in gene name', domain_name)
        return []
        # raise re.error
        
    intervals = get_intervals(idxs, window)
    for (start, end) in intervals:
        selected_text.append(all_text[start:end])
    return selected_text


def make_snippets(query, max_pages, snippet_window_size):
    sess = requests.Session()
    sess.mount('https://', HTTPAdapter(max_retries=retries))
    
    found_mentions = dict()
    try:
        # can return None
        data = get_first_page(query, sess)
    except requests.exceptions.RetryError:
        return None
    
    if data is None:
        return None
        
    if 'resultList' not in data:
        print(f'INFO: no data from EuropePMC for query {query}!')
        num_papers = -1
    else:
        num_papers = len(data['resultList']['result'])
    
    if num_papers <= 0:
        return found_mentions
    
    for k in range(max_pages):
        for i in range(len(data['resultList']['result'])):
            if 'fullTextUrlList' in data['resultList']['result'][i]:
                url_list =  data['resultList']['result'][i]['fullTextUrlList']['fullTextUrl']
                for elem in url_list:
                    doc_style = elem['documentStyle']
                    url = elem['url']
                    if doc_style == 'html' and url.find('europepmc.org') != -1:
                        paper_id = url.split('/')[-1]
                        try:
                            r = sess.get(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{paper_id}/fullTextXML')
                        except requests.exceptions.RetryError:
                            return None
                            
                        # should i work with other types of status?
                        if r.status_code == 200:
                            mentions = parse_xml_response(query, r.text, snippet_window_size)
                            found_mentions[paper_id] = mentions
                        
        if 'nextPageUrl' in data:
            try:
                rr = sess.get(data['nextPageUrl'])
            except requests.exceptions.RetryError:
                return None
            data = rr.json()
        else:
            break

    return found_mentions


def save_gene_stats(query, dir, genes):
    for protein_type in ['unreviewed', 'reviewed']:
        save_path = f'{dir}/{query}/genes_full_list_{protein_type}.txt'
        with open(save_path, 'w') as f:
            for gene_name in genes[protein_type]:
                print(gene_name, file=f)


def read_gene_list(path):
    gene_names = []
    with open(path, 'r') as f:
        for line in f:
            gene_name, uniprot_acc = line.strip().split('\t')
            gene_names.append((gene_name, uniprot_acc))
    return gene_names


def pull_genes(family, dir, max_pages, force_flag):
    print('LOG: Started pulling genes')
    # force_flag=True means the results will be re-recorded:
    # all gene names will be searched
    
    errors = 0
    database_name = get_dbname_from_acc(family)
    genes = dict({'reviewed': set(), 'unreviewed': set()})

    sess = requests.Session()
    sess.mount('https://', HTTPAdapter(max_retries=retries))

    for protein_type in ['unreviewed', 'reviewed']:
        print('LOG: Working with protein type:', protein_type)
        
        info_path = f'{dir}/{family}/genes_full_list_{protein_type}.txt'
        
        # if we don't want to rewrite results
        if (os.path.exists(info_path) and os.path.getsize(info_path) != 0 and not force_flag):
            print('INFO: path exists, not empty:', info_path)
            print('INFO: will be reading genes from this file. To rewrite genes, use flag FORCE=True.')
            genes[protein_type] = read_gene_list(info_path)
            continue
            
        url = f'https://www.ebi.ac.uk/interpro/api/protein/{protein_type}/entry/{database_name}/{family}'
        print('LOG: Starting retrieval')
        
        r = sess.get(url)
        
        if r.status_code == 204:
            print('INFO: Unsucessful: no data for protein type (error 204)', protein_type)
            continue
            
        if r.status_code != 200:
            print('ERROR: initial request unsuccessful! Stopping with this protein type:', protein_type)
            print('status code:', r.status_code)
            print(url)
            errors += 1
            continue
            
        data = json.loads(r.text)
        if 'results' not in data:
            print('ERROR: No key "results" in data! Stopping')
            errors += 1
            continue
            
        proteins = data['results']
        next_page_url = None
        if 'next' in data:
            next_page_url = data['next']
            
        for k in range(max_pages):
            for elem in proteins:
                gene_name = elem['metadata']['gene']
                uniprot_acc = elem['metadata']['accession']
                if gene_name is not None:
                    genes[protein_type].add((gene_name, uniprot_acc))
                    
            if next_page_url is None:
                break
                
            r = sess.get(next_page_url)
            if r.status_code != 200:
                print('ERROR: error while iterating through pages, stopping.')
                print('Status code:', r.status_code)
                print(next_page_url)
                errors += 1
                break
                
            data = json.loads(r.text)
            
            if 'results' not in data:
                print('ERROR: retrieved data has no key "results", stopping')
                errors += 1
                break
                
            proteins = data['results']
            
            if 'next' in data:
                next_page_url = data['next']
            else:
                next_page_url = None
                
        # maybe lowercase better?
        genes[protein_type] = set([(gene[0].casefold(), gene[1]) for gene in genes[protein_type]])
        # why similar genes appear? sometimes in reviewed and unreviewed are the same gene names

    for protein_type in ['unreviewed', 'reviewed']:
        info_path = f'{dir}/{family}/genes_full_list_{protein_type}.txt'
        
        if (force_flag or not os.path.exists(info_path) or os.path.getsize(info_path) == 0):
            print('LOG: writing genes_log to:', info_path)
            with open(info_path, 'w') as f:
                for gene in genes[protein_type]:
                    print(f'{gene[0]}\t{gene[1]}', file=f)
            print(f'INFO: wrote {len(genes[protein_type])} genes')

    if errors == 0:
        print('INFO: pull_genes function completed corectly')
    else:
        print('ERROR: pull_genes did not run correctly')
        
    return genes


# this is a supplementary function
# it was used to retrieve  uniprot_accs only for reviewed genes
def pull_reviewed_proteins(family, dir, max_pages=100):
    errors = 0
    database_name = get_dbname_from_acc(family)

    uniprot_accs = []

    sess = requests.Session()
    sess.mount('https://', HTTPAdapter(max_retries=retries))

    protein_type = 'reviewed'
    print('LOG: Working with type:', protein_type)
        
    url = f'https://www.ebi.ac.uk/interpro/api/protein/{protein_type}/entry/{database_name}/{family}'
    print('LOG: Starting retrieval')
    
    r = sess.get(url)
    
    if r.status_code == 204:
        print('INFO: Unsucessful: no data for protein type (error 204)', protein_type)
        return []
        
    if r.status_code != 200:
        print('ERROR: initial request unsuccessful! Stopping with this protein type:', protein_type)
        print('status code:', r.status_code)
        print(url)
        errors += 1
        return []
        
    data = json.loads(r.text)
    if 'results' not in data:
        print('ERROR: No key "results" in data! Stopping')
        errors += 1
        return []
        
    proteins = data['results']
    next_page_url = None
    if 'next' in data:
        next_page_url = data['next']
        
    for k in range(max_pages):
        for elem in proteins:
            uniprot_acc = elem['metadata']['accession']
            uniprot_accs.append(uniprot_acc)
                
        if next_page_url is None:
            break
            
        r = sess.get(next_page_url)
        if r.status_code != 200:
            print('ERROR: error while iterating through pages, stopping.')
            print('Status code:', r.status_code)
            print(next_page_url)
            errors += 1
            break
            
        data = json.loads(r.text)
        
        if 'results' not in data:
            print('ERROR: retrieved data has no key "results", stopping')
            errors += 1
            break
            
        proteins = data['results']
        
        if 'next' in data:
            next_page_url = data['next']
        else:
            next_page_url = None

        # if len(uniprot_accs) >= 100:
        #     break
                
        
    if errors == 0:
        print('INFO: pull_reviewed_proteins function completed corectly')
    else:
        print('ERROR: pull_reviewed_proteins did not run correctly')
        
    return uniprot_accs


def get_genes_by_type(query, dir_name, protein_type):
    result = set()
    genes_path = f'{dir_name}/{query}/genes_full_list_{protein_type}.txt'
    if not os.path.exists(genes_path):
        return result
    else:
        with open(genes_path, 'r') as f:
            for i, line in enumerate(f):
                try:
                    gene_name, uniprot_acc = line.strip().split('\t')
                except:
                    print(f'ERROR: error in processing {genes_path}. Cannot parse line {i}')
                result.add(gene_name)
    return result


def get_per_gene_snippets_summary(query, dir_name):
    data = dict()
    if os.path.exists(f'{dir_name}/{query}/genes_with_papers_list_all_review_status.txt'):
        filename = f'{dir_name}/{query}/genes_with_papers_list_all_review_status.txt'
        print('INFO: found all_review file', filename)
    elif os.path.exists(f'{dir_name}/{query}/genes_with_papers_list_unknown_review_status.txt'):
        filename = f'{dir_name}/{query}/genes_with_papers_list_unknown_review_status.txt'
        print('INFO: found unknown_review file', filename)
    else:
        raise Exception
        
    reviewed_genes = get_genes_by_type(query, dir_name, 'reviewed')
    unreviewed_genes = get_genes_by_type(query, dir_name, 'unreviewed')
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if len(line) != 3:
                continue
            gene_name, num_papers, num_snippets = line
            num_papers = int(num_papers)
            num_snippets = int(num_snippets)

            review_status = '0_unknown_review_status'
            if gene_name in unreviewed_genes:
                review_status = '1_unreviewed'
            if gene_name in reviewed_genes:
                review_status = '2_reviewed' 
                
            data[gene_name] = (num_papers, num_snippets, review_status)
            
    return data
    

def get_save_gene_snippets(query, dir_name, max_pages_per_gene, snippet_window_size, force_flag, from_gene_list=False, gene_list_filename=None, max_genes_each_type=3000):
    print('LOG: Started fetching gene snippets')
    
    log_file = f'{dir_name}/{query}/gene_snippet_search_log.json'
    snippets_path = f'{dir_name}/{query}/snippets_per_gene'

    def get_genes_from_txt(from_gene_list=False):
        genes = dict({'reviewed': set(), 'unreviewed': set(), 'unknown_review_status': set()})
    
        # if the data about genes already exists
        if from_gene_list:
            if gene_list_filename is None:
                print('ERROR: Gene list not provided. Either provide gene list or set from_gene_list=False.')
                raise Exception
                
            info_path = f'{dir_name}/{query}/{gene_list_filename}'
            with open(info_path, 'r') as f:
               for line in f:
                    gene_name, uniprot_acc = line.strip().split('\t')
                    genes['unknown_review_status'].add((gene_name.casefold(), uniprot_acc))
        else:
            for protein_type in ['reviewed', 'unreviewed']:
                info_path = f'{dir_name}/{query}/genes_full_list_{protein_type}.txt'
                with open(info_path, 'r') as f:
                    for line in f:
                        gene_name, uniprot_acc = line.strip().split('\t')
                        genes[protein_type].add((gene_name, uniprot_acc))
        return genes
        
    
    def log_progress(gene_name, status):
        try:
            with open(log_file, 'r') as f:
                log_data = json.load(f)
        except FileNotFoundError:
            log_data = {}
        
        log_data[gene_name] = status
        
        with open(log_file, 'w') as f:
            json.dump(log_data, f, indent=4)

    
    def verbose_progress():
        try:
            with open(log_file, 'r') as f:
                log_data = json.load(f)
                num_processed = len(log_data)
                num_success = list(log_data.values()).count('success')
                num_fail = list(log_data.values()).count('fail')
                print(f'INFO: tried to process: {num_processed}, success: {num_success}, fail: {num_fail}.')
        except:
            print('INFO: no gene snippet search log file yet')

    
    def get_processed_genes():
        try:
            with open(log_file, 'r') as f:
                log_data = json.load(f)
                return {gene_name for gene_name, status in log_data.items() if status == "success"}
        except FileNotFoundError:
            return set()

    def process_genes(unique_gene_names_all, max_pages_per_gene, snippet_window_size, max_genes_each_type):
        unique_gene_names_all = list(unique_gene_names_all)[:max_genes_each_type]
        # filter names that will cause error!
        # x.find('.') == -1 and
        unique_gene_names = [x for x in unique_gene_names_all if (x.find('/') == -1)]
        print(f'INFO: removed {len(unique_gene_names_all) - len(unique_gene_names)} genes with forward slash (/) in name.')
        
        processed_genes = get_processed_genes()
        already_processed = 0
    
        for gene_name in unique_gene_names:
            output_path = snippets_path + '/' + gene_name + '.csv'
            if os.path.exists(output_path):
                # need to update status? I made a separate script for this
                # print(f'INFO: Skipping {gene_name}: snippets already exist.')
                already_processed += 1
                continue
            if gene_name in processed_genes:
                # print(f"INFO: Skipping {gene_name}, already processed. Empty.")
                already_processed += 1
                continue
            
            found_mentions = make_snippets(query=gene_name,
                                           max_pages=max_pages_per_gene,
                                           snippet_window_size=snippet_window_size)
            
            if found_mentions is not None:
                log_progress(gene_name, "success")
            else:
                log_progress(gene_name, "fail")
                continue
            
            num_snippets = sum([len(found_mentions[paper_id]) for paper_id in found_mentions])
            if num_snippets > 0:
                # can cause error btw :)
                make_and_save_data(gene_name, output_path=output_path, found_mentions=found_mentions)
            else:
                print(f'INFO: for {gene_name} found 0 snippets')
        print(f'INFO: genes that were already processed for snippets: {already_processed}')

    
    def get_gene_snippet_info(gene_name):
        num_papers, num_snippets = 0, 0
        df_path = f'{dir_name}/{query}/snippets_per_gene/{gene_name}.csv'
        if os.path.exists(df_path):
            df = pd.read_csv(df_path)
            num_papers = len(set(df['paper_id']))
            num_snippets = len(df)
        return num_papers, num_snippets

    
    def make_summary():
        summary_info = []
        gene_names_all = [x[:-4] for x in os.listdir(snippets_path)]
        for gene_name in gene_names_all:
            num_papers, num_snippets = get_gene_snippet_info(gene_name)
            if num_papers > 0 and num_snippets > 0:
                summary_info.append((gene_name, num_papers, num_snippets))
            
        with open(f'{dir_name}/{query}/genes_with_papers_list_all_review_status.txt', 'w') as f:
            for elem in summary_info:
                print(f'{elem[0]}\t{elem[1]}\t{elem[2]}', file=f)
        print('INFO: created gene snippet search summary in format: gene_name\tnum_papers\tnum_snippets')

    
    genes = get_genes_from_txt(from_gene_list)
    
    if (os.path.exists(snippets_path) and len(os.listdir(snippets_path)) != 0):
        if not force_flag:
            num_genes_with_snippets = len(os.listdir(snippets_path))
            print(f'Snippets already exist ({num_genes_with_snippets} files). Not making snippets again')
            return num_genes_with_snippets
        else:
            print('INFO: Snippets exist, but forced to make them again. Info: new files can be created, existing files will not be altered.')

    mkdirsafe(snippets_path)
    verbose_progress()
    
    for protein_type in ['reviewed', 'unreviewed', 'unknown_review_status']:
        print(f"INFO: Working with: {protein_type} proteins")
        
        if len(genes[protein_type]) != 0:
            unique_gene_names = sorted(list(set([x[0] for x in genes[protein_type]])))
            random.seed(24)
            random.shuffle(unique_gene_names)
            
            print(f'INFO: found {len(unique_gene_names)} unique gene names of type {protein_type}')
            process_genes(unique_gene_names, max_pages_per_gene, snippet_window_size, max_genes_each_type)
    
    print('INFO: gene snippet search completed successfully')

    make_summary()
    
    num_genes_with_snippets = len(os.listdir(snippets_path))
    print(f'INFO: Total {len(os.listdir(snippets_path))} genes with snippets.')
    
    return num_genes_with_snippets


def enumerate_snippets(query, dir_name):
    # renaming columns
    def adjust_column_names(df):
        try:
            df.rename(columns={'excerpt_num': 'snippet_num'}, inplace=True)
            df.drop(columns=['excerpt_id'], inplace=True)
        except:
            pass
        try:
            df.drop(columns=['snippet_id'], inplace=True)
        except:
            pass
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
            
        return df
    
    snippets_path = f'{dir_name}/{query}/snippets_per_gene'
    file_names = os.listdir(snippets_path)
    file_names = [x for x in file_names if x[0] != '.']

    paper_to_cnt = defaultdict(int)
    snippet_id_to_gene = dict()
    paper_to_genes = dict()

    for file_name in file_names:
        snippet_file_path = snippets_path + '/' + file_name
        df = pd.read_csv(snippet_file_path)
        df = adjust_column_names(df)
        
        gene_name = os.path.splitext(file_name)[0]
        snippet_ids = []
        
        for i, row in df.iterrows():
            paper_id = row['paper_id']
            
            if paper_id not in paper_to_genes:
                paper_to_genes[paper_id] = []
            paper_to_genes[paper_id].append(gene_name)
            
            paper_snippet_counter = paper_to_cnt[paper_id]
            paper_to_cnt[paper_id] += 1
            
            snippet_id = f'{paper_id}_{paper_snippet_counter}'
            snippet_ids.append(snippet_id)
            
            snippet_id_to_gene[snippet_id] = gene_name
            
        df.insert(2, 'snippet_id', snippet_ids)
        df.to_csv(snippet_file_path)
            
    with open(f"{dir_name}/{query}/tmp/snippet_id_to_gene.json", "w") as json_file:
        json.dump(snippet_id_to_gene, json_file, indent=4)
    with open(f"{dir_name}/{query}/tmp/paper_to_genes.json", "w") as json_file:
        json.dump(paper_to_genes, json_file, indent=4)
    paper_to_cnt = dict(paper_to_cnt)
    with open(f"{dir_name}/{query}/tmp/paper_to_cnt.json", "w") as json_file:
        json.dump(paper_to_cnt, json_file, indent=4)

    print('LOG: Sucessfully added snippet IDs. snippet_id_to_gene.josn, paper_to_genes.json, paper_to_cnt.json created')


def get_snippet_id_to_gene_mapping(query, dir_name):
    with open(f"{dir_name}/{query}/tmp/snippet_id_to_gene.json", "r") as json_file:
        snippet_id_to_gene = json.load(json_file)
    return snippet_id_to_gene
    

def run_factcheck(client, text1, text2, config):
    pass
    
    
def remove_citations(text):
    pattern = r'\s*\[\d+(?:,\s*\d+)*(?:-\d+)?\]'
    cleaned_text = re.sub(pattern, '', text)
    cleaned_text = re.sub(r'\s+', ' ', cleaned_text).strip()
    cleaned_text = cleaned_text.replace('\n', '')
    return cleaned_text


def join_snippets_into_prompt(query, dir_name, run_name, N, config):
    snippets_path = f'{dir_name}/{query}/snippets_per_gene'
    save_path = f'{dir_name}/{query}/{run_name}/per_gene_joined_prompts_dirty'
        
    mkdirsafe(save_path)
    
    file_names = os.listdir(snippets_path)
    file_names = [x for x in file_names if x[0] != '.']
    
    for file_name in file_names:
        gene_name = os.path.splitext(file_name)[0]
        df = pd.read_csv(snippets_path + '/' + file_name)
        
        instr_1 = config['instr_1_gene_temp'].format(gene_name)
        
        if config['instr_2_gene_temp'].find('file:') != -1:
            with open(config['instr_2_gene_temp'][5:], 'r') as f:
                instr_2_fam_temp = f.read()
        else:
            instr_2_fam_temp = config['instr_2_gene_temp']
            
        instr_2 = instr_2_fam_temp.format(gene_name)
    
        prompts, n_chars, n_snippets = [], [], []
        part = ''
        i = 0

        for j, row in df.iterrows():
            snippet, snippet_id = row['snippet'], row['snippet_id']
            snippet = remove_citations(snippet)
            i += 1
            data = '{' + 'ID: ' + snippet_id + ',\n' + 'Content: ' + snippet + '}\n'
            part += data
            
            if i == N:
                prompt = instr_1 + '\n' + part + '\n' + instr_2
                prompts.append(prompt)
                n_chars.append(len(prompt))
                n_snippets.append(N)
                part = ''
                i = 0
        
        if len(part) > 0:
            prompt = instr_1 + '\n' + part  + '\n' + instr_2
            prompts.append(prompt)
            n_chars.append(len(prompt))
            n_snippets.append(i)
        
        df_prompts = pd.DataFrame({'prompt': prompts,
                                   'num_chars': n_chars,
                                   'num_snippets': n_snippets
                                  })

        
        fname = f'{save_path}/{gene_name}.csv'
        df_prompts.to_csv(fname)

    print('INFO: Joined snippets into prompts using config')


def select_genes(query, dir_name, run_name, N, TH_GOOD=0.5, TH_BAD=0.15, prompt_th=0, fam_filter_th=0):
    def fam_filter(gene_names, snippet_stats):
        result = []
        gene_stat_dict = dict()
        
        selected_genes_filename = f'{dir_name}/{query}/{run_name}/selected_genes.txt'
        selected_genes_file = open(selected_genes_filename, 'w')
    
        spec_stats_filename = f'{dir_name}/{query}/tmp/gene_spec_stats.tsv'
        # if not os.path.exists(spec_stats_filename):
        make_spec_stats_file(dir_name, query)
        
        with open(spec_stats_filename, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                gene_name, all_hits_total, value_true, value_false = line
                all_hits_total = int(all_hits_total)
                value_true = float(value_true)
                value_false = float(value_false)
                gene_stat_dict[gene_name] = (all_hits_total, value_true, value_false)
            
        for gene_name in gene_names:
            if gene_name not in gene_stat_dict:
                continue
            all_hits_total, value_true, value_false = gene_stat_dict[gene_name]
            
            if (value_true >= 0.5 and value_false <= 0.2) or (value_false < 0.01 and value_true > 0.1):
                result.append(gene_name)
                print(f'{gene_name}\t{all_hits_total}\t{value_true}\t{value_false}\t{snippet_stats[gene_name][0]}\t{snippet_stats[gene_name][1]}\t{snippet_stats[gene_name][2]}', file=selected_genes_file)
                if len(result) >= N:
                    break
                
        selected_genes_file.close()
        print('INFO: Created gene selection file', selected_genes_filename)
        return result

    snippet_stats = get_per_gene_snippets_summary(query, dir_name)
    print(f'INFO: {len(snippet_stats)} genes with snippets were found. Doing specificity filtering...')
    if len(snippet_stats) == 0:
        return []
    
    gene_names, gene_info = zip(*snippet_stats.items())
    num_snippets = [x[1] for x in gene_info]
    gene_type = [x[2] for x in gene_info]
    _, __, gene_names = zip(*sorted(zip(gene_type, num_snippets, gene_names), reverse=True))
    
    gene_names = fam_filter(gene_names, snippet_stats)
    gene_names = gene_names[:N]

    print(f'INFO: Selected {len(gene_names)} genes:', gene_names)
    file_paths = [x + '.csv' for x in gene_names]
    
    return file_paths


def get_selected_genes_filepaths(query, dir_name, run_name, N):
    gene_names = []
    selected_genes_filename = f'{dir_name}/{query}/{run_name}/selected_genes.txt'
    with open(selected_genes_filename, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gene_name = line[0]
            gene_names.append(gene_name)

    file_paths = [x + '.csv' for x in gene_names]
    file_paths = file_paths[:N]
    
    return file_paths
    

def find_citations(text):
    ''' 
    Finds all citations in the text.
    '''
    pattern = r'\[\w+_\w+\]'
    matches = re.findall(pattern, text)
    
    matches = [x[1:-1] for x in matches]
    return matches


def get_snippet_dict_by_gene_name(query, dir_name, gene_name):
    ''' Returns dictionary snippet_id -> snippet text for a specified gene name'''
    
    snippets_path = f'{dir_name}/{query}/snippets_per_gene/{gene_name}.csv'
    if not os.path.exists(snippets_path):
        raise Exception(f"Snippet path error: file {snippets_path} doesnt exist")
        
    df_snippets = pd.read_csv(snippets_path)
    snippet_dict = dict(zip(df_snippets['snippet_id'], df_snippets['snippet']))
    return snippet_dict


def make_factcheck_df(query, dir_name, snippet_ids_, gene_names_, snippet_content_, type='unk'):
    df = pd.DataFrame({'snippet_id': snippet_ids_,
                        'gene_name': gene_names_,
                        'snippet': snippet_content_,
                        })
    if len(df) == 0:
        print('WARNING: no citations found. Not creating factchech dataframe. This should not happen normally!')
        return df
    
    df['paper_id'] = df['snippet_id'].str.split('_').str[0]
    df.insert(1, 'paper_id', df.pop('paper_id'))
    df['snippet_error'] = df['snippet'].str.startswith('ERROR', na=False)
    
    df.drop_duplicates(inplace=True)
    df.to_csv(f'{dir_name}/{query}/factcheck_{type}.csv')
    return df
    

def get_gene_summary_citations(query, dir_name, response_df_path):
    # todo: counts for each paper
    ''' 
    Reads per-gene summaries, finds citations in them,
    prints citation id + snippet in the output csv file.
    Requires summary path to exist!
    '''
    snippet_ids_, gene_names_, snippet_content_ = [], [], []
    
    df_summaries = pd.read_csv(response_df_path)
    for i, row in df_summaries.iterrows():
        gene_name, text = row['gene_name'], row['response']
        
        snippet_dict = get_snippet_dict_by_gene_name(query, dir_name, gene_name)
        citations = find_citations(text)

        for snippet_id in citations:
            snippet_ids_.append(snippet_id)
            gene_names_.append(gene_name)
            if snippet_id in snippet_dict:
                snippet_content_.append(snippet_dict[snippet_id])
            else:
                snippet_content_.append('ERROR: False snippet ID!')

    df_factcheck = make_factcheck_df(query, dir_name, snippet_ids_, gene_names_, snippet_content_, type='genes')
    return df_factcheck


def get_family_summary_citations(query, dir_name, text):
    '''
    Requires summary path to exist!
    '''
    snippet_ids_, gene_names_, snippet_content_ = [], [], []
    snippet_id_to_gene = get_snippet_id_to_gene_mapping(query, dir_name)
    
    citations = find_citations(text)
    
    for snippet_id in citations:
        snippet_ids_.append(snippet_id)
        
        if snippet_id not in snippet_id_to_gene:
            gene_names_.append('ERROR: hallucinated snippet ID')
            snippet_content_.append('ERROR: hallucinated snippet ID')
        else:
            gene_name = snippet_id_to_gene[snippet_id]
            gene_names_.append(gene_name)
            snippet_dict = get_snippet_dict_by_gene_name(query, dir_name, gene_name)
            snippet_content_.append(snippet_dict[snippet_id])
            
    df_factcheck = make_factcheck_df(query, dir_name, snippet_ids_, gene_names_, snippet_content_, type='family')
    if len(df_factcheck) == 0:
        return df_factcheck
        
    print('INFO: Logging used snippets per gene for this family summary:')
    print(df_factcheck['gene_name'].value_counts())

    error_df = df_factcheck[df_factcheck['snippet'].str.startswith('ERROR', na=False)]
    error_count = error_df.shape[0]
    print(f'INFO: detected {error_count} bad (hallucinated) snippet ids.')
    if error_count > 0:
        print('INFO: Erros are below')
        print(error_df)

    return df_factcheck
    

def copy_response(query, dir_name, save_path):
    pass
