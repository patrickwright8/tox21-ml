import pandas as pd
import numpy as np
import os

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

# Setup I/O
root = os.getcwd() 
test_datapath = os.path.join(root, 'projects\\tox21-ml\\data\\test_data')

eval_results_path = os.path.join(test_datapath, 'tox21_eval_dataset_results.csv')
eval_smiles_path = os.path.join(test_datapath, 'tox21_eval_dataset_smiles.csv')
test_sdf_path = os.path.join(test_datapath, 'tox21_test_dataset_results.txt')

train_set_path = os.path.join(root, 'projects\\tox21-ml\\data\\combined_tox21.csv')

export_path = os.path.join(test_datapath, 'combined_test_set.csv')

# Defines RDKit fragment cleaner
largest_Fragment = rdMolStandardize.LargestFragmentChooser()
# %% Generate dataframe with test smiles (cleaned) and assay results

# Read (and clean) smiles and assay results from SDF file
test_smiles = []
test_results = []
data = Chem.SDMolSupplier(test_sdf_path)

for mol in data:
    if mol == None:
        continue
    assay_results = mol.GetPropsAsDict()
    test_results.append(assay_results)
    largest_mol = largest_Fragment.choose(mol) # Smiles are cleaned here
    sm = Chem.MolToSmiles(largest_mol) 
    test_smiles.append(sm)

# Convert assay results to dataframe form
test_results_df = pd.DataFrame(test_results).iloc[:,2:]

# Convert smiles from list form to pd series
test_smiles_df = pd.DataFrame(test_smiles, columns=['smiles'])

# Combines smiles with assay results; drops duplicates and resets index
test_dfs_uncombined = [test_smiles_df, test_results_df]
test_set = pd.concat(test_dfs_uncombined, axis=1).drop_duplicates(subset=['smiles']).reset_index(drop=True)

# %% Generate dataframe with eval smiles and assay results

# Read eval smiles from CSV
eval_smiles = pd.read_csv(eval_smiles_path, names=['smiles', 'Sample ID'], delimiter='\t', header=0)

eval_results = pd.read_csv(eval_results_path, delimiter='\t').replace('x', np.nan)

eval_smiles_list = list(eval_smiles.iloc[:,0])
eval_smiles_cleaned = []

for sm in eval_smiles_list:
    m = Chem.MolFromSmiles(sm)
    if m != None:
        largest_mol = largest_Fragment.choose(m)
        sm = Chem.MolToSmiles(largest_mol)
        eval_smiles_cleaned.append(sm)
    if m == None:
        eval_smiles_cleaned.append(np.nan)

eval_smiles.iloc[:,0] = eval_smiles_cleaned
eval_smiles = eval_smiles.dropna(inplace=False).drop_duplicates(inplace=False).reset_index(drop=True)

eval_set = pd.merge(eval_smiles, eval_results, on='Sample ID')
eval_set = eval_set.drop(columns='Sample ID')
# %% Combine dataframes into one test set! Export as CSV

combined_test_set = pd.concat([test_set, eval_set]).reset_index(drop=True)
combined_test_set = combined_test_set.drop_duplicates(subset=['smiles'])

# Generates set of all smiles present in training set
train_set = set(pd.read_csv(train_set_path).iloc[:,0])

# Removes rows in the test set containing smiles present in the training set
test_set_filtered = combined_test_set[~combined_test_set['smiles'].isin(train_set)]

# Shuffles test set and resets index
final_test_set = test_set_filtered.sample(frac=1).reset_index(drop=True)

# Exports to CSV
final_test_set.to_csv(export_path, index=False)
