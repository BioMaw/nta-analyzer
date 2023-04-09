#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   main.py
@Time    :   2023/04/09 10:49:57
@Author  :   Bio-Maw 
@Version :   1.0
@Contact :   martint.maw@gmail.com
@Desc    :   None
'''

import pandas as pd
from pathlib import Path
import seaborn as sns

def parse_file(file):
    #read contents of tsv
    df = pd.read_csv(file)  

convert_to_float = lambda x: float(x.replace(',', '.'))


INPUT_DIR = Path('/mnt/c/Projects/Data/ev_sizes')
FIRST_LINE =74
# generate folder for headerless files
NO_HEAD_DIR = Path(INPUT_DIR, 'no_head')
NO_HEAD_DIR.mkdir(exist_ok=True)

file_list= list(INPUT_DIR.glob('*.txt'))
df_master = pd.DataFrame()

for file in file_list:
    # remove header from file and put it in headerless folder
    fn_no_head = Path(NO_HEAD_DIR, file.name) 
    with open(file, "r") as fin:
        rows = fin.readlines()[FIRST_LINE:]
    with open(fn_no_head, 'w') as fout:
        fout.writelines(rows)

    #read headerless files
    df = pd.read_csv(fn_no_head, sep ='\t') 
    # convert to numeric values 
    for column in df.columns:
        df[column] = df[column].apply( lambda x: float(x.replace(',', '.')))
        df[column] = df[column].astype(float)
    # Get indices of rows to drop (-1 row) and drop them
    drop_indices = df[df['Number'] == -1].index
    df = df.drop(index=drop_indices)

    # read information out of filename (be careful to name correctly!!)
    part_list = str(fn_no_head.stem).split('_')
    df['biol_replic']=part_list[3]
    df['treatment']=part_list[4]
    df['dilution']=int(part_list[5][3:])
    df['tech_replic']=int(part_list[8])

    # calculate additional columns
    df['conc_corr'] = df['Concentration / cm-3'] * df['dilution']
    # drop redundant columns
    df = df.drop(columns=['Number', 'Volume / nm^3','Area / nm^2','Concentration / cm-3'])

    # append to master dataframe
    df_master = pd.concat([df_master, df])
df_master.reset_index(drop=True, inplace=True)

#clean up some errors in User Input, make all treatment lower case
df_master['treatment']=df_master['treatment'].str.lower()

# normalize each technical replicate by their own total concentration
sum_conc = df_master.groupby(['biol_replic', 'tech_replic', 'treatment'])['conc_corr'].transform('sum')

# normalization to 1
df_master['conc_corr_norm'] = df_master['conc_corr'] / sum_conc

## averaging of technical replic
df_avg = df_master.groupby(['biol_replic', 'treatment', 'Size / nm']).mean()
df_avg.reset_index(inplace=True)
df_avg.drop(columns='tech_replic', inplace=True)

means = df_avg.groupby(['biol_replic', 'treatment']).mean()
## violinplot improvement needed
sns.violinplot(data=df_avg, x= 'biol_replic', y='conc_corr_norm', hue ='treatment', split=True)