import os
import json
import pandas as pd
import pickle
from utils import mkdirsafe


def get_num_genes(dir_name, acc):
    genes_file = f'{dir_name}/{acc}/gene_list_sorted.txt'
    
    genes_set = set()
    with open(genes_file, 'r') as f:
        for line in f:
            gene_name, uniprot_acc = line.strip().split('\t')
            genes_set.add(gene_name.casefold())
    n = len(genes_set)
    n_bad = sum([1 for x in genes_set if (x.find('/') != -1)])
    return n, n_bad
    

def get_num_genes_with_snippets(dir_name, acc):
    path = f'{dir_name}/{acc}/snippets_per_gene'
    if not os.path.exists(path):
        print('INFO: snippets dir does not exist')
        return 0
    l = len(os.listdir(path))
    return l


def get_log_stats(dir_name, acc):
    log_file = f'{dir_name}/{acc}/gene_snippet_search_log.json'
    with open(log_file, 'r') as f:
        log_data = json.load(f)
        num_processed = len(log_data)
        num_success = list(log_data.values()).count('success')
        num_fail = list(log_data.values()).count('fail')
    return num_processed, num_success, num_fail


def verbose_search_stats(dir_name, acc):
    num_processed, num_success, num_fail = get_log_stats(dir_name, acc)
    n_genes, n_bad_genes = get_num_genes(dir_name, acc)
    n_genes_with_snippets = get_num_snippets(dir_name, acc)
    
    print(acc, 'Genes:', n_genes, 'Bad (unsuitable names) genes:', n_bad_genes, 'Processed:', num_processed, 'With snippets:', n_genes_with_snippets)


def checker_for_rerun(dir_name, acc, verbose=False):
    num_processed, num_success, num_fail = get_log_stats(acc)
    n_genes, n_bad_genes = get_num_genes(acc)
    n_good_genes = n_genes - n_bad_genes
    n_genes_with_snippets = get_num_snippets(acc)
    
    if (n_good_genes < 1000 and num_processed < n_good_genes) or \
    (n_good_genes >= 1000 and num_processed < 1000 - n_bad_genes):
        if verbose:
            print('Rerunning:', acc, n_good_genes, num_processed, '---', n_genes_with_snippets)
        return True
    return False


def checker_too_few_snippets(dir_name, acc, num_genes_max,cmin_genes_with_snippets=10, verbose=False):
    n_genes, n_bad_genes = get_num_genes(acc)
    n_genes_with_snippets = get_num_snippets(acc)
    if n_genes_with_snippets < min_genes_with_snippets and n_genes - n_bad_genes > num_genes_max - 1000:
        if verbose:
            print(acc, 'with snippets:', n_genes_with_snippets, 'total good genes:', n_genes - n_bad_genes)
        return True
    return False