# %% Imports
import chemprop

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG
    
# %% Hyperparameter optimization

arguments = [
    '--data_path', 'C:/projects/tox21-ml/data/combined_tox21.csv',
    '--dataset_type', 'classification',
    '--save_dir', 'C:/projects/tox21-ml/ckpts',
    '--gpu', '0',
    '--batch_size', '50',
    '--smiles_columns', 'Smiles',
    '--target_columns', 'nr-ahr', 'nr-ar-lbd', 'nr-ar', 'nr-aromatase',\
        'nr-er-lbd', 'nr-er', 'nr-ppar-gamma', 'sr-are', 'sr-atad5', 'sr-hse',\
            'sr-mmp', 'sr-p53',]

extra_args = ['--num_iters', '100', 
    '--search_parameter_keywords', 'basic',
    '--config_save_path', 'C:/projects/tox21-ml/ckpts',
    '--manual_trial_dirs', 'C:/projects/tox21-ml/old_ckpts/trial_seed_0',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_1',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_2',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_3',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_4',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_5',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_6',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_7',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_8',\
    'C:/projects/tox21-ml/old_ckpts/trial_seed_9',\
    '--startup_random_iters', '10'
]

### Notes:
    # Suppose 1 iteration (30 epochs) takes ~30 minutes. Then for
    # just 20 iterations, bayesian optimization would take 
    # ~10 hours on this dataset! 
    # Best to run overnight or over the weekend.
    # Still - worth it to squeeze out every last bit of performance.

hyperopt_args = arguments + extra_args

hyperopt_args = chemprop.args.HyperoptArgs().parse_args(hyperopt_args)
chemprop.hyperparameter_optimization.hyperopt(hyperopt_args)

# %% Train the model using the winning parameters

# Coming soon!