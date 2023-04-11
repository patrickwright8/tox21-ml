# %% Imports
import chemprop
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG


root = os.path.join(os.getcwd(), 'projects\\tox21-ml')
datapath = os.path.join(root, 'data\\combined_tox21.csv')
savedir = os.path.join(root, 'train_ckpts')
configpath = os.path.join(root, 'hyperopt\\best.json')
    
# %% Perform 10-fold CV on single model

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
    '--num_folds', '10',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--split_type', 'random',
    '--config_path', configpath,
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',
]

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)

# %% Run 10-fold CV on an ensemble of 10 models using the winning hyperparameters

arguments = [
    '--data_path', datapath,
    '--dataset_type', 'classification',
    '--save_dir', savedir,
    '--gpu', '0',
    '--batch_size', '50',
    '--num_folds', '10',
    '--ensemble_size', '10'
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--split_type', 'random',
    '--config_path', configpath,
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',
]

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)