# import csv
import os
# import numpy as np
import pandas as pd
import config as cfg
# import json
# import matplotlib.pyplot as plt
# import seaborn as sns


def merge_terms_features():
    df_features = pd.read_csv(cfg.data['mobidb'] + 'mobidb_regions_residues_filtered.csv', sep='\t')
    df_terms = pd.read_csv(cfg.data['disprot'] + 'disprot_regions_residues_expanded.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)},)

    features = ['curated-disorder-disprot', 'prediction-low_complexity-merge', 'prediction-lip-anchor', 'prediction-proline_rich-mobidb_lite_sub',
                     'prediction-negative_polyelectrolyte-mobidb_lite_sub', 'prediction-positive_polyelectrolyte-mobidb_lite_sub',
                     'prediction-polyampholyte-mobidb_lite_sub', 'prediction-polar-mobidb_lite_sub']
    df_merge = pd.merge(df_terms, df_features, how='outer', on=['UniProt_id', 'residue'])
    #
    for feature in features:
        count = 0
        binary_list = []
        for residue in df_merge.iterrows():
            print(count)
            count += 1
            if residue[1].feature == feature:
                binary_list.append(1)
            else:
                binary_list.append(0)
        df_merge[feature] = binary_list

    df_merge = df_merge.drop(columns= 'feature')
    #
    file = cfg.data['merged'] + 'terms_features.csv'
    df_merge.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)



    print(df_terms)

def expanded_sequence_length():
    df = pd.read_csv(cfg.data['merged'] + 'terms_features.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})
    df_seq = pd.read_csv(cfg.data['merged'] + 'terms_features.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})[['UniProt_id','sequence_length']].drop_duplicates()
    df_seq = df_seq.loc[(df_seq['sequence_length'] != '')]

    df.drop('sequence_length', 1, inplace=True)

    df_merge = pd.merge(df, df_seq, how='inner', on=['UniProt_id'])
    file = cfg.data['merged'] + 'terms_features_length.csv'
    df_merge.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)

def compress_features():
    df = pd.read_csv(cfg.data['merged'] + 'terms_features_length.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})
    df = df.fillna(-1)
    df.loc[df['term_id'] == '', 'term_id'] = -1
    aggregation_functions = { 'UniProt_id': 'first', 'DisProt_id': 'first',  'residue': 'first', 'term_id': 'first', 'sequence_length': 'first', 'curated-disorder-disprot': 'sum', 'prediction-low_complexity-merge': 'sum', 'prediction-lip-anchor': 'sum', 'prediction-proline_rich-mobidb_lite_sub': 'sum', 'prediction-negative_polyelectrolyte-mobidb_lite_sub': 'sum', 'prediction-positive_polyelectrolyte-mobidb_lite_sub': 'sum', 'prediction-negative_polyelectrolyte-mobidb_lite_sub': 'sum', 'prediction-polyampholyte-mobidb_lite_sub': 'sum', 'prediction-polar-mobidb_lite_sub': 'sum'}
    df_new = df.groupby(['UniProt_id', 'residue', 'DisProt_id', 'term_id']).aggregate(aggregation_functions)
    file = cfg.data['merged'] + 'terms_features_compressed.csv'
    df_new.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)


def add_sequence():
    df = pd.read_csv(cfg.data['merged'] + 'terms_features_compressed.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})
    df_sequences = pd.read_csv(cfg.data['disprot'] + 'disprot_residues.csv', sep='\t',  converters={'sequence_length': lambda x: str(x)})
    df_merge = pd.merge(df, df_sequences, how='inner', on=['UniProt_id', 'residue'])

    file = cfg.data['merged'] + 'terms_features_w_amino.csv'
    df_merge.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)


def add_group():
    df = pd.read_csv(cfg.data['merged'] + 'terms_features_w_amino.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})
    df_groups = pd.read_csv(cfg.data['functions'] + 'groups.csv', sep='\t',  converters={'IDPO_term_id': lambda x: str(x), 'GO_term_id': lambda x: str(x)})
    df_merge = pd.merge(df, df_groups[['term_id','group']], how='inner', on=['term_id'])
    df_unique = df_merge.drop_duplicates(subset=['UniProt_id' , 'residue', 'group'])

    group_list = df_unique['group'].drop_duplicates()

    total_amino_counts = pd.DataFrame()
    for group in group_list:

        df_group = df_unique.loc[(df_unique['group'] == group)]
        amino_count = df_group.groupby(["UniProt_id", "char"]).apply(count_mean).reset_index(name='size')
        amino_count['group'] = group
        total_amino_counts = total_amino_counts.append(amino_count, ignore_index=True)
    file = cfg.data['functions'] + 'amino_groups.csv'
    total_amino_counts.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)


def count_mean(group):
    for i, row in group.iterrows():
        return group.shape[0]
