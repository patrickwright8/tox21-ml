# %% Imports
import chemprop
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = os.getcwd()
datapath = os.path.join(root, 'data\\combined_tox21.csv')
test_path = os.path.join(root, 'data\\test_data\\eval_set.csv')
savedir = os.path.join(root, 'hyperopt_ckpts')
hyperparam_savedir = os.path.join(savedir, 'best.json')
    
# %% Hyperparameter optimization

# NOTE:
    # For scafold split to work, you must change "np.float" in 
    # C:\ProgramData\Anaconda3\envs\chemprop\Lib\site-packages\chemprop\data
    # to "float" or "np.float64"

arguments = [
    '--data_path', datapath,
    '--dataset_type', 'classification',
    '--save_dir', savedir,
    '--gpu', '0',
    '--batch_size', '50',
    '--num_folds', '3',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--split_type', 'scaffold_balanced',
    '--split_size', '.9', '.1', '0',
    '--separate_test_path', test_path,
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',]

extra_args = ['--num_iters', '100', 
    '--search_parameter_keywords', 'basic',
    '--config_save_path', hyperparam_savedir,
]

hyperopt_args = arguments + extra_args

hyperopt_args = chemprop.args.HyperoptArgs().parse_args(hyperopt_args)
chemprop.hyperparameter_optimization.hyperopt(hyperopt_args)