# %% Imports, I/O
import chemprop
import pandas as pd
import numpy as np 
import pickle as pkl
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = 'C:/projects/tox21-ml' # Eventually change to os.getcwd()
datafolder = root + '/data/'
datapath = datafolder + 'combined_tox21.pkl'

with open(datapath, 'rb') as r:
    df = pkl.load(r)
    
# %% Hyperparameter optimization

arguments = [
    '--data_path', 'C:/projects/chemprop/datasets/train.csv',
    '--dataset_type', 'regression',
    '--save_dir', 'ckpts',
    '--gpu', '0',
    '--batch_size', '4',
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',]

extra_args = ['--num_iters', '20', # Only for hyperopt
    '--search_parameter_keywords', 'basic', # Only for hyperopt
    '--config_save_path', 'C:/projects/chemprop/tox21-ckpts'
]

hyperopt_args = arguments + extra_args

hyperopt_args = chemprop.args.HyperoptArgs().parse_args(hyperopt_args)
chemprop.hyperparameter_optimization.hyperopt(hyperopt_args)