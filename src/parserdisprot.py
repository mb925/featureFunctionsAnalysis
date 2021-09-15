import os
import pandas as pd
import config as cfg
import json
import numpy as np
from pandas import read_csv

def filter_tsv():
    print(cfg.data['disprot'])
    df = read_csv(cfg.data['disprot'] + 'disprot_2021.csv', sep='\t')[['acc', 'disprot_id', 'start', 'end', 'term', 'ec', 'region_sequence']]
    df['sequence_length'] = df['region_sequence'].apply(lambda x: len(x))
    df = df.drop(['region_sequence'], axis=1)
    df['term_id'] = df['term'].apply(lambda x: x.split(':')[2])
    df['ec_id'] = df['ec'].apply(lambda x: x.split(':')[2])
    df = df.drop(['term'], axis=1)
    df = df.drop(['ec'], axis=1)
    df = df.drop_duplicates(subset=['acc', 'disprot_id', 'start', 'end', 'term_id']) # in case there are duplicates because of different pubmed id
    df.rename(columns={'acc': 'UniProt_id', 'disprot_id': 'DisProt_id'}, inplace=True)
    file = cfg.data['disprot'] + 'disprot_regions.csv'
    df.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)

    print(df)

def json_to_csv():
    with open(cfg.data['disprot'] + 'disprot_regions.json') as json_file:
        data = json.load(json_file)
        df = pd.DataFrame()
        for el in data['data']:

            sequence_length = el['length']
            uniprot_id = el['acc']
            disprot_id = el['disprot_id']

            for region in el['regions']:
                if region['term_namespace'] == 'Disorder function':
                    df_el = pd.DataFrame({
                        'UniProt_id': uniprot_id,
                        'DisProt_id': disprot_id,
                        'term_id': str(region['term_id']),
                        'ec_id': [str(region['ec_id'])],
                        'start': [region['start']],
                        'end': [region['end']],
                        'sequence_length': [sequence_length]
                    })
                    df = df.append(df_el)
                    print(df)

        file = cfg.data['disprot'] + 'disprot_regions.csv'
        df.to_csv(file, mode='a',  header=(not os.path.exists(file)), sep='\t', index=False)

def all_uniprot_residues():
    df = read_csv(cfg.data['disprot'] + 'disprot_regions.csv', sep='\t', converters={'term_id': lambda x: str(x)})
    uniprots = {}
    uniprots_all_residues = pd.DataFrame({'UniProt_id': [], 'residue': []})
    uniprots_id = []
    residues = []
    for row in df.iterrows():
        uniprot_id = row[1]['UniProt_id']
        uniprots[uniprot_id] = row[1]['sequence_length']
    for uniprot, length in uniprots.items():
        for residue in range(1, length + 1):
            uniprots_id.append(uniprot)
            residues.append(residue)

    uniprots_all_residues['UniProt_id'] = uniprots_id
    uniprots_all_residues['residue'] = residues

    df_regions = read_csv(cfg.data['disprot'] + 'disprot_regions_residues.csv', sep='\t', converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})

    df_expanded_regions = pd.merge(uniprots_all_residues, df_regions, how='left', on=['UniProt_id', 'residue']).drop_duplicates()
    file = cfg.data['disprot'] + 'disprot_regions_residues_expanded.csv'
    df_expanded_regions.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)
    print(uniprots)




def start_end_to_residue():
    df = read_csv(cfg.data['disprot'] + 'disprot_regions.csv', converters={'term_id': lambda x: str(x), 'ec_id': lambda x: str(x)}, sep='\t')
    df_residues = pd.DataFrame()

    for i, row in df.iterrows():

        for r in range(row.start, row.end + 1):
            df_el = pd.DataFrame({
                'UniProt_id': row['UniProt_id'],
                'DisProt_id': row['DisProt_id'],
                'term_id': row['term_id'],
                'residue': [r],
                'sequence_length': [row['sequence_length']]
            })
            df_residues = df_residues.append(df_el)
            print(df_el)
    file = cfg.data['disprot'] + 'disprot_regions_residues.csv'
    df_residues.to_csv(file, mode='a',  header=(not os.path.exists(file)), sep='\t', index=False)

