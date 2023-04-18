# %% Imports, I/O
import chemprop
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = os.path.join(os.getcwd(), 'projects\\tox21-ml')
datapath = os.path.join(root, 'data\\test_data\\combined_test_set.csv')
# configpath = os.path.join(root, 'hyperopt\\best.json')
ensemble_savedir = os.path.join(root, 'ensemble_ckpts')

# %% Run ensemble predictions on test set



# %% Evaluates predictions by comparing to true values

