import os
import pandas as pd
import config as cfg
from pandas import read_csv



list_features = ['curated-disorder-disprot','prediction-low_complexity-merge', 'prediction-lip-anchor', 'prediction-proline_rich-mobidb_lite_sub',
                     'prediction-negative_polyelectrolyte-mobidb_lite_sub', 'prediction-positive_polyelectrolyte-mobidb_lite_sub',
                     'prediction-polyampholyte-mobidb_lite_sub', 'prediction-polar-mobidb_lite_sub']


def create_file_residues():
    df = pd.read_csv(cfg.data['mobidb'] + 'mobidb_regions.tsv', sep='\t')
    df_disprot = pd.read_csv(cfg.data['disprot'] + 'disprot_regions.csv', sep='\t')
    # df = df.loc[df.feature == 'prediction-low_complexity-merge']
    L = df_disprot['UniProt_id']
    df = df[~df.acc.isin(L)]

    df_residues = pd.DataFrame()

    for i, row in df.iterrows():
        print(row)
        start_ends = row['start..end'].split(',')
        for region in start_ends:
            start = int(region.split('..')[0])
            end = int(region.split('..')[1])
            for residue in range(start, end + 1):
                df_el = pd.DataFrame({
                    'UniProt_id': row['acc'],
                    'feature': row['feature'],
                    'residue': [residue],
                    'sequence_length': row['length']
                })
                df_residues = df_residues.append(df_el)
                print(df_el)
    file = cfg.data['mobidb'] + 'mobidb_regions_residues.csv'
    df_residues.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)


def create_dicts():

    df_length = pd.read_csv(cfg.data['mobidb'] + 'mobidb_regions.tsv', sep='\t')
    df_mobi = pd.read_csv(cfg.data['mobidb'] + 'mobidb_regions_residues.csv', sep='\t')

    dicts = {}
    for feature in list_features:

        df_mobi_ft = df_mobi.loc[df_mobi.feature == feature]

    uniprots = set(df_mobi_ft['UniProt_id'])
    for uniprot in uniprots:
        len = df_length.loc[df_length['UniProt_id'] == uniprot].sequence_length
        for res in range(1, len):
            if df_mobi_ft.loc[df_mobi_ft['UniProt_id'] == uniprot & df_mobi_ft.residue == res].empty:
                dicts[feature][uniprot + '_' + res] = 0
            else:
                dicts[feature][uniprot + '_' + res] = 1

    file = cfg.data['mobidb'] + 'mobidb_features.csv'
    # df_features.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)


def filter_residues_by_disprot_id():
    uniprot_list = read_csv(cfg.data['disprot'] + 'disprot_regions.csv', sep='\t', converters={'term_id': lambda x: str(x)})['UniProt_id'].drop_duplicates()
    df_mobi = pd.read_csv(cfg.data['mobidb'] + 'mobidb_regions_residues.csv', sep='\t')


    boolean_series = df_mobi.UniProt_id.isin(uniprot_list)
    filtered_df = df_mobi[boolean_series]
    file = cfg.data['mobidb'] + 'mobidb_regions_residues_filtered.csv'
    filtered_df.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)
    print(filtered_df)



