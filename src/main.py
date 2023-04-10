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
import matplotlib.pyplot as plt

def parse_file(file):
    #read headerless files
    df = pd.read_csv(file, sep ='\t') 
    # convert to numeric values 
    for column in df.columns:
        df[column] = df[column].apply( lambda x: float(x.replace(',', '.')))
        df[column] = df[column].astype(float)
    # Get indices of rows to drop (-1 row) and drop them
    drop_indices = df[df['Number'] == -1].index
    df = df.drop(index=drop_indices)

    # read information out of filename (be careful to name correctly!!)
    part_list = str(file.stem).split('_')
    df['biol_replic']=part_list[3]
    df['treatment']=part_list[4]
    df['dilution']=int(part_list[5][3:])
    df['tech_replic']=int(part_list[8])

    # calculate additional columns
    df['conc_corr'] = df['Concentration / cm-3'] * df['dilution']
    # drop redundant columns
    df = df.drop(columns=['Number', 'Volume / nm^3','Area / nm^2','Concentration / cm-3'])
    return df
INPUT_DIR = Path('/mnt/c/Projects/Data/ev_sizes')
FIRST_LINE =74
# generate folder for headerless files
NO_HEAD_DIR = Path(INPUT_DIR, 'no_head')
NO_HEAD_DIR.mkdir(exist_ok=True)

file_list= list(INPUT_DIR.glob('*.txt'))
df_master = pd.DataFrame()

for file in file_list:
    # remove header from file and put it in headerless folder
    #TODO check if headerless already exists
    fn_no_head = Path(NO_HEAD_DIR, file.name) 
    with open(file, "r") as fin:
        rows = fin.readlines()[FIRST_LINE:]
    with open(fn_no_head, 'w') as fout:
        fout.writelines(rows)
    df = parse_file(fn_no_head)
    # append to master dataframe
    df_master = pd.concat([df_master, df])
df_master.reset_index(drop=True, inplace=True)

#clean up some errors in User Input, make all treatment lower case
df_master['treatment']=df_master['treatment'].str.lower()

### start dataset manipulation
# truncate based on size
# cut off above the highest occuring size
#max_valid_size = max(df_master.loc[df_master['conc_corr'] != 0, 'Size / nm'].unique())
max_valid_size = 800
df_master=df_master.loc[df_master['Size / nm'] < max_valid_size]

# cut off below 30 nm 
min_cutoff = 30
df_master=df_master.loc[df_master['Size / nm'] > min_cutoff]

# add classification based on size
is_small_class = (df_master['Size / nm'] > 30) & (df_master['Size / nm'] <= 150) 
is_large_class = (df_master['Size / nm'] >= 151) & (df_master['Size / nm'] <= 800)
df_master.loc[is_small_class, 'Size_Class'] = 'small_class'
df_master.loc[is_large_class, 'Size_Class'] = 'large_class'

# normalize each technical replicate by their own total concentration
# get sum for each treatment of each tech_replic of each biol_replic
sum_conc = df_master.groupby(['biol_replic', 'tech_replic', 'treatment'])['conc_corr'].transform('sum')

# normalization to 1
df_master['conc_corr_norm'] = df_master['conc_corr'] / sum_conc

## averaging over technical replic
df_avg = df_master.groupby(['biol_replic', 'treatment', 'Size / nm']).mean()
df_avg.reset_index(inplace=True)
df_avg.drop(columns='tech_replic', inplace=True)

# Calculate mean within each group
mean_df = df_avg.groupby(['biol_replic', 'treatment'])[['Size / nm', 'conc_corr_norm']].apply(
    lambda x: (x['Size / nm'] * x['conc_corr_norm']).sum() / x['conc_corr_norm'].sum())

# Rename columns
mean_df = mean_df.reset_index()
mean_df.columns = ['biol_replic', 'treatment', 'Mean']

# plot histograms:
# Get list of unique cell lines
cell_lines = df_avg['biol_replic'].unique()

# Create subplot grid
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))

for i, cell_line in enumerate(cell_lines):
    row = i // 2
    col = i % 2
    ax = axes[row, col]

    p_data = df_avg.loc[df_avg['biol_replic'] == cell_line]
    sns.histplot(data=p_data, x='Size / nm', weights='conc_corr_norm', hue='treatment', multiple='stack', binwidth=5, ax=ax)
    ax.set_title('Cell line {}'.format(cell_line))

# Add overall title and axis labels
plt.suptitle('Histograms by cell line and treatment')
fig.text(0.5, 0.04, 'Size / nm', ha='center')
fig.text(0.04, 0.5, 'Count', va='center', rotation='vertical')

# Adjust spacing between plots
plt.tight_layout()

# Show plot
plt.show()


# means = df_avg.groupby(['biol_replic', 'treatment']).mean()
# mean_value = (df['Size / nm'] * df['conc_corr_norm']).sum() / df['Count'].sum()
# ## violinplot improvement needed
# sns.violinplot(data=df_avg, x= 'biol_replic', y='conc_corr_norm', hue ='treatment', split=True)