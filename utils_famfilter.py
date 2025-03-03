import requests
import os

def get_cross_references(family, data):
    ''' Works with one page of the results '''
    
    all_hits, false_hits, pos_hits = 0, 0, 0
    all_hits = len(data['results'])
    
    for item in data['results']:
        has_correct, has_incorrect = 0, 0
        
        for info in item['uniProtKBCrossReferences']:
            db_name, id_, props = info['database'], info['id'], info['properties']
            if db_name == 'PANTHER':
                # all_hits += 1
                if id_.find(':') != -1:
                    continue
                # id_ = id_.split(':')[0]

                if id_ == family:
                    has_correct += 1
                else:
                    has_incorrect += 1
        if has_correct > 0:
            pos_hits += 1
        if (has_incorrect > 0 and has_correct == 0):
            false_hits += 1
    
    return all_hits, pos_hits, false_hits


def get_hits_for_gene_name(family, gene_name, num_pages=20, verb=0):
    all_proteins = 0
    all_hits_total, pos_hits_total, false_hits_total  = 0, 0, 0
    
    query = f"https://rest.uniprot.org/uniprotkb/search?&query=gene:{gene_name}&size=100"
    
    for i in range(num_pages):        
        r = requests.get(query)
        data = r.json()
        all_proteins += len(data['results'])
        
        all_hits, pos_hits, false_hits = get_cross_references(family, data)
        all_hits_total += all_hits
        false_hits_total += false_hits
        pos_hits_total += pos_hits
        
        if 'Link' in r.headers:
            link_text = r.headers['Link']
            link_text = link_text.split(';')[0]
            query = link_text[1:-1]
        else:
            break

    return all_hits_total, pos_hits_total, false_hits_total


def get_stats_for_gene_name(family, gene_name, verb_errors=False):
    try:
        all_hits_total, pos_hits_total, false_hits_total  = get_hits_for_gene_name(family, gene_name)
    except:
        if verb_errors:
            print(f'INFO: Error in retrieving data for gene name {gene_name}')
        raise requests.exceptions.RequestException
        
    if all_hits_total == 0:
        if verb_errors:
            print(f'INFO: no information about genes for family {family} retrieved.')
        raise ZeroDivisionError
    
    prop_true = pos_hits_total/all_hits_total
    prop_false = false_hits_total/all_hits_total
    
    return all_hits_total, prop_true, prop_false


def make_spec_stats_file(dir_name, family):
    # from scratch
    
    summary_path = f'{dir_name}/{family}/tmp/gene_spec_stats.tsv'
    print('INFO: Creating gene specificity file', summary_path)
    
    summ_file = open(summary_path, 'w')

    gene_names = os.listdir(f'{dir_name}/{family}/snippets_per_gene')
    gene_names = [x[:-4] for x in gene_names]
    print(f'INFO: Found {len(gene_names)} gene names')
    
    for gene_name in gene_names:
        try:
            all_hits_total, prop_true, prop_false = get_stats_for_gene_name(family, gene_name)
        except:
            print('INFO: Error in retrieving data for gene name', gene_name)
            continue
            
        if all_hits_total == 0:
            print('INFO:', family, 'Division by zero attempted')
            continue

        # sort first?
        print(f'{gene_name}\t{all_hits_total}\t{round(prop_true, 2)}\t{round(prop_false, 2)}', file=summ_file)
        print(f'{gene_name}\t{all_hits_total}\t{round(prop_true, 2)}\t{round(prop_false, 2)}')

    summ_file.close()
