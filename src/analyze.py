# import csv
import os
# import numpy as np
import pandas as pd
import config as cfg

def intesection_count(group):
    count = 0
    for i, row in group.iterrows():
        if row['term_id'] == '00012' and row['prediction-low_complexity-merge'] == 1:
            count +=1

    return count


def union_count(group):
    count = 0
    for i, row in group.iterrows():
        if row['term_id'] == '00012' or row['prediction-low_complexity-merge'] == 1:
            count += 1

    return count


def intersection_union():
    # test - example: prediction-low_complexity-merge, 00012
    df = pd.read_csv(cfg.data['merged'] + 'terms_features_length.csv', sep='\t',  converters={'term_id': lambda x: str(x), 'sequence_length': lambda x: str(x)})
    df = df.loc[(df['prediction-low_complexity-merge'] == 1) | (df['term_id'] == '00012')]
    res_int = df.groupby('UniProt_id').apply(intesection_count)
    res_int.name = 'intersection'
    res_int = res_int.reset_index()
    res_uni = df.groupby('UniProt_id').apply(union_count)
    res_uni.name = 'union'
    res_uni = res_uni.reset_index()
    df_merge = pd.merge(res_int, res_uni, how='inner', on=['UniProt_id'])
    file = cfg.data['analyze'] + 'union_intersection.csv'
    df_merge.to_csv(file, mode='a', header=(not os.path.exists(file)), sep='\t', index=False)
    print(df)


