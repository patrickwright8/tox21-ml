import pandas as pd
import numpy as np
import os
import pickle as pkl

root = 'C:/projects/tox21-ml'
datafolder = root + '/data/'
filenames = os.listdir(datafolder)
filenames = [i for i in filenames if i.count('.smiles')>0]

dfs = []
smiles = np.nan
for file in filenames:
    temp_filepath = datafolder + file
    assay_name = file.replace('.smiles','')
    temp_df = pd.read_csv(temp_filepath, sep='\s+|\t', engine='python',
                          header=None,
                          names = ['Smiles', 'ID', str(assay_name)],
                          usecols = [0,2])
    if smiles == np.nan:
        smiles = temp_df.iloc[:,0]
    dfs.append(temp_df)
    
def merge_and_drop_duplicates(df_list):
    merged = pd.DataFrame()
    for df in df_list:
        if merged.empty:
            merged = df
        else:
            merged = pd.merge(merged, df, on='Smiles', how='outer')
            merged = merged.drop_duplicates()
    return merged

# dfs = dfs[:2]

# Test!
df = merge_and_drop_duplicates(dfs)

savename = datafolder + 'combined_tox21.pkl'
with open(savename, 'wb') as w:
    pkl.dump(df, w)

# %% PD tests

import pandas as pd


# create some sample dataframes
df1 = pd.DataFrame({'smiles': ['C1=CC=CC=C1', 'CC(=O)O'], 'value1': [1, 2]})
df2 = pd.DataFrame({'smiles': ['C1=CC=CC=C1', 'CCN'], 'value2': [3, 4]})
df3 = pd.DataFrame({'smiles': ['CCN', 'CC(=O)O'], 'value3': [5, 6]})

merged_df = pd.DataFrame({'smiles': []})  # create an empty dataframe with only the 'smiles' column

dfs = [df1, df2, df3]

for df in dfs:
    merged_df = pd.merge(merged_df, df, on='smiles', how='outer')  


# %% 
# merge dataframes based on 'smiles' column
merged_df = df1.merge(df2, on='smiles', how='outer') \
              .merge(df3, on='smiles', how='outer')

print(merged_df)









