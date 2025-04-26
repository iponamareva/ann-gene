''' This is a script that receives a family accession,
    scrapes the list of genes belonging to this family,
    tries to find papers mentioning each gene,
    summarize gene function first and then
    make a summary of the family.

    This script uses OpenAI GPT-4 API,
    which means it will cost some money to run it.
    Be careful with number of families and parameters
    you run it with. 
'''


import os
import json
import argparse
import pandas as pd
import openai

from datetime import datetime

from utils import verbose_args, get_first_page, mkdirsafe, make_and_save_data, make_snippets, pull_genes, get_save_gene_snippets, join_snippets_into_prompt, copy_response
from utils import save_args_log, enumerate_snippets
from utils import select_genes
from utils import pull_genes_for_uniprot

from utils_gpt import get_gpt_genes_response, get_gpt_family_response, verbose_gpu_usage
from utils_gpt import factcheck_summary, factcheck_gene_summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', '--mode', type=str, required=True, choices=["from-fam-acc", "from-gene-list", "from-uniprot-list"], default="from-fam-acc")
    
    # parser.add_argument('-config', '--config', type=str, default='configs/config-desc-name-CoT.json')
    parser.add_argument("-q", "--query", type=str, help='family name; actual name for gene search or placeholder name, depending on the mode')
    parser.add_argument('-dir', '--dir-name', type=str, default='output_per_query')
    parser.add_argument('-o', '--text-output-dir-name', type=str, default='output')

    parser.add_argument('-uniprot-list', '--uniprot-list', type=str, default="uniprot_list.txt")
    parser.add_argument("---gene-list", '--gene-list', type=str, default='gene_list_sorted.txt')
    
    parser.add_argument("-F", "--FORCE", type=bool, default=False)
    
    parser.add_argument("-max", "--max-pages-per-family", type=int, default=1000, help='maximum number of pages per family to parse, default=1000')
    parser.add_argument("-max2", "--max-pages-per-gene", type=int, default=1, help='maximum number of pages per gene to parse, default=1')
    parser.add_argument("-max3", "--max-genes-each-type", type=int, default=1000, help='maximum number of genes')
    parser.add_argument('-s', '--snippet-window-size', type=int, default=300, help='snippet window size in characters, snippet will be 2 times longer')
    
    parser.add_argument('-N1', '--num-snippets-in-prompt', type=int, default=100)
    parser.add_argument('-N', '--gpt4-n', type=int, default=10)
    parser.add_argument('-run-gpt', '--run-gpt', type=int, default=0)
    parser.add_argument('-sf', '--do-spec-filter', type=bool, default=True)
    
    parser.add_argument('-v', '--verbose', type=int, default=0)

    args = parser.parse_args()
    verbose_args(args)

    # with open(args.config, 'r') as config_file:
    #     config = json.load(config_file)
    with open("/hps/software/users/agb/research/irina/configs/config-desc-name-CoT.json", 'r') as config_file:
        config = json.load(config_file)

    mkdirsafe(args.text_output_dir_name)
    # raise Exception(f'ERROR: Text piyput dir {args.text_output_dir_name} does not exist. Please create it')

    mkdirsafe(f'{args.dir_name}')
    mkdirsafe(f'{args.dir_name}/{args.query}')
    mkdirsafe(f'{args.dir_name}/{args.query}/tmp')
    now = datetime.now()
    run_name = now.strftime(f"run_results_{config['model']}_%Y-%m-%d %H:%M:%S")
    mkdirsafe(f'{args.dir_name}/{args.query}/{run_name}')
    save_args_log(args, config, run_name)
    print(f'LOG: run_name {run_name}. Arguments logged in run_args.log')

    ''' Getting and saving snippets; Joining snippets into prompts '''

    if args.mode == "from-fam-acc":
        pull_genes(family=args.query, dir=args.dir_name, max_pages=args.max_pages_per_family, force_flag=args.FORCE)
        num_genes_with_snippets = get_save_gene_snippets(args.query, args.dir_name, args.max_pages_per_gene, args.snippet_window_size, force_flag=args.FORCE, from_gene_list=False, max_genes_each_type=args.max_genes_each_type, gene_list_filename=args.gene_list)

    if args.mode == "from-gene-list":
        num_genes_with_snippets = get_save_gene_snippets(args.query, args.dir_name, args.max_pages_per_gene, args.snippet_window_size, force_flag=args.FORCE, from_gene_list=True, max_genes_each_type=args.max_genes_each_type, gene_list_filename=args.gene_list)
    
    if args.mode == "from-uniprot-list":
        gene_list_filename = pull_genes_for_uniprot(family=args.query, dir_name=args.dir_name, uniprot_accs_path=args.uniprot_list, max_pages=args.max_pages_per_family, force_flag=args.FORCE)
        num_genes_with_snippets = get_save_gene_snippets(args.query, args.dir_name, args.max_pages_per_gene, args.snippet_window_size, force_flag=args.FORCE, from_gene_list=True, max_genes_each_type=args.max_genes_each_type, gene_list_filename=gene_list_filename)

    enumerate_snippets(args.query, args.dir_name)
    
    join_snippets_into_prompt(args.query, args.dir_name, run_name, args.num_snippets_in_prompt, config)
    selected_genes_file_paths = select_genes(args.query, args.dir_name, run_name, args.gpt4_n, spec_filter=args.do_spec_filter)
    
    ''' Using GPT-4 API to make summaries '''
    if args.run_gpt:
        if len(selected_genes_file_paths) > 0:
    
            client = openai.OpenAI()
            GPT_USAGE_1, gene_names, gpt4_responses = get_gpt_genes_response(client, args.query, args.dir_name, run_name, args.gpt4_n, config)
            GPT_USAGE_2, parsed_response, pmid_parsed_response = get_gpt_family_response(client, args.query, args.dir_name, run_name, gene_names, gpt4_responses, config, args.text_output_dir_name)
    
            print('*'*30 + args.query + '*'*30, f'\n{parsed_response}\n', '*'*79)
            print(f'\n{pmid_parsed_response}\n', '*'*79)
        else:
            print('INFO: No genes with snippets were selected. Family summary is not possible')

if __name__ == "__main__":
    main()