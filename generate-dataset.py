import pandas as pd
import numpy as np
import os
import pickle as pkl

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem

root = 'C:/projects/tox21-ml' # Eventually change to os.getcwd()
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

df = merge_and_drop_duplicates(dfs) # Merges assay results; drops duplicates
df = df.reset_index(drop=True) # Resets index

# Still need to clean smiles by removing salts, extra molecule fragments
smiles = list(df.iloc[:,0])
largest_Fragment = rdMolStandardize.LargestFragmentChooser()
for count, sm in enumerate(smiles):
    try:
        m = Chem.MolFromSmiles(sm)
    except: 
        pass
    if m == None:
        smiles[count] = np.nan
        continue
    largest_mol = largest_Fragment.choose(m)
    new_sm = Chem.MolToSmiles(largest_mol)
    smiles[count] = new_sm
smiles = pd.Series(smiles, name='Smiles')
df.iloc[:,0] = smiles
df = df.dropna(subset='Smiles').reset_index(drop=True)

savename = datafolder + 'combined_tox21.pkl'
with open(savename, 'wb') as w:
    pkl.dump(df, w)

print('\nSaved cleaned, merged DataFrame to file!\n')