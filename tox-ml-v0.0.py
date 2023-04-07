# %% Imports
import chemprop

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG
    
# %% Hyperparameter optimization

# NOTE:
    # For scafold split to work, you must change "np.float" in 
    # C:\ProgramData\Anaconda3\envs\chemprop\Lib\site-packages\chemprop\data
    # to "float" or "np.float64"

arguments = [
    '--data_path', 'C:/Users/ptw80/projects/tox21-ml/data/combined_tox21.csv',
    '--dataset_type', 'classification',
    '--save_dir', 'C:/Users/ptw80/projects/tox21-ml/ckpts',
    '--gpu', '0',
    '--batch_size', '50',
    '--num_folds', '3',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--split_type', 'scaffold_balanced',
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',]

extra_args = ['--num_iters', '500', 
    '--search_parameter_keywords', 'basic',
    '--config_save_path', 'C:/Users/ptw80/projects/tox21-ml/ckpts/best.json',
]

hyperopt_args = arguments + extra_args

hyperopt_args = chemprop.args.HyperoptArgs().parse_args(hyperopt_args)
chemprop.hyperparameter_optimization.hyperopt(hyperopt_args)

# %% Train the model using the winning parameters

arguments = [
    '--data_path', 'C:/Users/ptw80/projects/tox21-ml/data/combined_tox21.csv',
    '--dataset_type', 'classification',
    '--save_dir', 'C:/Users/ptw80/projects/tox21-ml/ckpts/random_results',
    '--num_folds', '10',
    '--gpu', '0',
    '--batch_size', '50',
    '--depth', '3',
    '--dropout', '.3',
    '--ffn_num_layers', '3',
    '--hidden_size', '900',
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53']

args = chemprop.args.TrainArgs().parse_args(arguments)
mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)