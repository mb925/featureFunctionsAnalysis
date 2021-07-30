# import csv
import os
# import numpy as np
import pandas as pd
import config as cfg
import json
# import matplotlib.pyplot as plt
# import seaborn as sns

def json_to_tsv():
    with open(cfg.data['data'] + 'disprot_regions.json') as json_file:
        data = json.load(json_file)
        dict_terms = {}
        df = pd.DataFrame()

        for el in data['data']:

            sequence_length = el['length']
            uniprot_id = el['acc']
            disprot_id = el['disprot_id']

            for region in el['regions']:
                term_id = region['term_id']
                ec_id = region['ec_id']
                start = region['start']
                df_el = pd.DataFrame({
                    'UniProt_id': uniprot_id,
                    'DisProt_id': disprot_id,
                    'term_id': region['term_id'],
                    'ec_id': region['ec_id'],
                    'start': region['start'],
                    'end': region['end'],
                    'sequence_length': sequence_length
                })
                df = df.append(df_el)

            file = cfg.data['data'] + 'disprot_regions.csv'
            df.to_csv(file, mode='a',  header=(not os.path.exists(file)), sep='\t', index=False)
