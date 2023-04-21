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
import math
import openpyxl
import numpy as np


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
    df['ID']=int(part_list[2])
    df['biol_replic']=part_list[3]
    df['treatment']=part_list[4]
    df['dilution']=int(part_list[5][3:])
    #print('dilution for file' + str(file.stem) + ' = ' + str(int(part_list[5][3:]))) 
    df['tech_replic']=int(part_list[8])

    # calculate additional columns
    df['conc_corr'] = df['Concentration / cm-3'] * df['dilution']
    # drop redundant columns
    df = df.drop(columns=['Number', 'Volume / nm^3','Area / nm^2','Concentration / cm-3'])
    return df

def plot_bars(df, export = False):
    #plot histograms:
    # Get list of unique cell lines
    cell_lines = df['biol_replic'].unique()
    nrows = math.floor(math.sqrt(len(cell_lines)))
    ncols = math.ceil(math.sqrt(len(cell_lines)))
    # Create subplot grid
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 8),
                             sharey=True)
    
    for i, cell_line in enumerate(cell_lines):
        row = i // nrows
        col = i % ncols
        ax = axes[row, col]

        p_data = df.loc[df['biol_replic'] == cell_line]
        sns.barplot(data=p_data, x='treatment', y='EV_count/Cell_count', ax=ax)
        ax.set_title('Cell line {}'.format(cell_line))
    # Add overall title and axis labels
    figtitle='Barplots by cell line and treatment EVcount divided by Cellcount.'
    plt.suptitle(figtitle)

    # Adjust spacing between plots
    plt.tight_layout()

    # Show plot
    plt.show()
    if export:
       fig.savefig(Path(export, figtitle + ' .png'), format='png', dpi=300, bbox_inches='tight') 


def plot_hist(df, binwidth = 5, export = False, binexport = False):
    #plot histograms:
    # Get list of unique cell lines
    cell_lines = df['biol_replic'].unique()
    nrows = math.floor(math.sqrt(len(cell_lines)))
    ncols = math.ceil(math.sqrt(len(cell_lines)))
    # Create subplot grid
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 8))

    for i, cell_line in enumerate(cell_lines):
        row = i // nrows
        col = i % ncols
        ax = axes[row, col]

        p_data = df.loc[df['biol_replic'] == cell_line]
        sns.histplot(data=p_data, x='Size / nm', weights='conc_corr_norm', 
                    hue='treatment',  binwidth=binwidth, ax=ax)
        ax.set_title('Cell line {}'.format(cell_line))

    # Add overall title and axis labels
    figtitle='Histograms by cell line and treatment. Binsize = ' + str(binwidth) + ' nm'
    plt.suptitle(figtitle)
    fig.text(0.5, 0.04, 'Size / nm', ha='center')
    fig.text(0.04, 0.5, 'Probability Mass Function', va='center', rotation='vertical')

    # Adjust spacing between plots
    plt.tight_layout()

    # Show plot
    plt.show()
    if export:
       fig.savefig(Path(export, figtitle + ' .png'), format='png', dpi=300, bbox_inches='tight') 
    if binexport:
        # Compute the histogram using numpy
        bins = np.arange(p_data['Size / nm'].min(), p_data['Size / nm'].max() + binwidth, binwidth)
        counts, bin_edges = np.histogram(p_data['Size / nm'], bins=bins, weights=p_data['conc_corr_norm'])
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_magnitudes = counts * binwidth

        # Convert the bin centers and magnitudes to a pandas dataframe
        bin_data = pd.DataFrame({'bin_center': bin_centers, 'bin_magnitude': bin_magnitudes})
        bin_data.to_excel(Path(export, 'bin_data_Binsize' + str(binwidth) + 'nm.xlsx'))

INPUT_DIR = Path('/mnt/c/Projects/Data/230421_ev_sizes')
FIRST_LINE =74
# generate folder for headerless files
NO_HEAD_DIR = Path(INPUT_DIR, 'no_head')
NO_HEAD_DIR.mkdir(exist_ok=True)

OUTPUT_DIR = Path(INPUT_DIR, 'exports')
OUTPUT_DIR.mkdir(exist_ok=True)

file_list= list(INPUT_DIR.glob('*.txt'))
df_master = pd.DataFrame()

xlsx_list= list(INPUT_DIR.glob('*.xlsx'))

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

for file in xlsx_list:
    df_desc = pd.read_excel(file)

df_master = pd.merge(df_master, df_desc, on='ID', suffixes=('', '_desc'))
df_master = df_master.drop(
    columns=[col for col in df_master.columns if col.endswith('_desc')])

#clean up some errors in User Input, make all treatment lower case
df_master['treatment']=df_master['treatment'].str.lower()

# remove 200 dilution data
df_master = df_master.loc[df_master['dilution'] != 200]

### start dataset manipulation
# truncate based on size
# cut off above the highest occuring size
#max_valid_size = max(df_master.loc[df_master['conc_corr'] != 0, 'Size / nm'].unique())
max_valid_size = 800
df_master=df_master.loc[df_master['Size / nm'] < max_valid_size]

# cut off below 30 nm 
min_cutoff = 30
df_master=df_master.loc[df_master['Size / nm'] > min_cutoff]

# normalize each technical replicate by their own total concentration
# first step: get sum of EVs across all sizes
# for each treatment of each tech_replic of each biol_replic
sum_conc = df_master.groupby(['biol_replic', 'tech_replic', 'treatment'])['conc_corr'].transform('sum')

# second step: normalization to 1
df_master['sum'] = sum_conc
df_master['conc_corr_norm'] = df_master['conc_corr'] / sum_conc

### check if number of nonzero sizes is equal for all treatments and replics
def count_unique_nonzero(x):
    return x[x != 0].nunique()

# Apply custom function to grouped DataFrame
unique_sizes = df_master.groupby(['biol_replic', 'tech_replic', 'treatment'])['conc_corr'].apply(count_unique_nonzero)
# Reset index to turn resulting Series into a DataFrame
unique_concentrations = unique_sizes.reset_index()

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

# add classification based on size
is_small_class = (df_avg['Size / nm'] > 30) & (df_avg['Size / nm'] <= 150) 
is_large_class = (df_avg['Size / nm'] >= 151) & (df_avg['Size / nm'] <= 800)
df_avg.loc[is_small_class, 'Size_Class'] = 'small_class'
df_avg.loc[is_large_class, 'Size_Class'] = 'large_class'

# add normalized EV count by Cellcount 
df_count_EV = df_master.groupby(['biol_replic', 'tech_replic', 'treatment'])[
    'Living Events/μL(V)','sum'].mean()
df_count_EV = df_count_EV.reset_index()
df_count_EV['EV/Cellcount'] = df_count_EV['sum'] / df_count_EV['Living Events/μL(V)']
#get average and std
df_count_stat = pd.DataFrame()
df_count_stat['mean'] = df_count_EV.groupby(['biol_replic', 'treatment'])['EV/Cellcount'].mean()
df_count_stat['std'] = df_count_EV.groupby(['biol_replic', 'treatment'])['EV/Cellcount'].std()
df_count_stat = df_count_stat.rename(columns={'mean': 'EV_count/Cell_count'})
df_count_stat = df_count_stat.reset_index()

# sns.barplot(x='treatment', y='mean', hue='biol_replic', data=df_count_stat)
plot_bars(df_count_stat, export = OUTPUT_DIR)

df_count_stat.to_excel(Path(OUTPUT_DIR, 'data_EV_by_cell_mean_std.xlsx'))
df_count_EV.to_excel(Path(OUTPUT_DIR, 'data_EV_by_cell.xlsx'))

plot_hist(df_avg, binwidth=2.5, export = OUTPUT_DIR, binexport = True)
plot_hist(df_avg, binwidth=5, export = OUTPUT_DIR, binexport = True)
plot_hist(df_avg, binwidth=10, export = OUTPUT_DIR, binexport = True)
plot_hist(df_avg, binwidth=20, export = OUTPUT_DIR, binexport = True)
plot_hist(df_avg, binwidth=50, export = OUTPUT_DIR, binexport = True)
plot_hist(df_avg, binwidth=120, export = OUTPUT_DIR, binexport = True)

#export table to excel
df_avg.to_excel(Path(OUTPUT_DIR, 'data_long_table.xlsx'))

#export mean results to excel
mean_df.to_excel(Path(OUTPUT_DIR, 'mean_dist_result.xlsx'))
# means = df_avg.groupby(['biol_replic', 'treatment']).mean()
# mean_value = (df['Size / nm'] * df['conc_corr_norm']).sum() / df['Count'].sum()
# ## violinplot improvement needed
# sns.violinplot(data=df_avg, x= 'biol_replic', y='conc_corr_norm', hue ='treatment', split=True)