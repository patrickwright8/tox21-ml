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

export_path = os.path.join(test_datapath, 'combined_test_set.csv')

# %% Generate dataframe with test smiles (cleaned) and assay results

# Read (and clean) smiles and assay results from SDF file
test_smiles = []
test_results = []
data = Chem.SDMolSupplier(test_sdf_path)
largest_Fragment = rdMolStandardize.LargestFragmentChooser()
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




# %% Combine dataframes into one test set and clean!




# %% Export cleaned dataset


# Don't forget to shuffle!
