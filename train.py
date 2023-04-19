# %% Imports
import chemprop
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG


root = os.getcwd()
datapath = os.path.join(root, 'data\\combined_tox21.csv')
configpath = os.path.join(root, 'hyperopt_results\\best.json')
ensemble_savedir = os.path.join(root, 'ensemble_ckpts_with_testing')
test_set_path = os.path.join(root, 'data\\test_data\\formatted_combined_test_set.csv')
    
# %% Run 5-fold CV on an ensemble of 3 models using the winning hyperparameters

arguments = [
    '--data_path', datapath,
    '--dataset_type', 'classification',
    '--save_dir', ensemble_savedir,
    '--gpu', '0',
    '--batch_size', '50',
    '--num_folds', '10',
    '--ensemble_size', '10',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--split_type', 'random',
    '--split_size', '.9', '.1', '0',
    '--separate_test_set', test_set_path,
    '--config_path', configpath,
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',
]

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)