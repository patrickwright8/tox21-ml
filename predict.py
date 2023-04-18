# %% Imports, I/O
import chemprop
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = os.path.join(os.getcwd(), 'projects\\tox21-ml')
datapath = os.path.join(root, 'data\\test_data\\combined_test_set.csv')
ensemble_savedir = os.path.join(root, 'ensemble_ckpts')
prediction_path = os.path.join(root, 'predict\\predictions.csv')

# %% Run ensemble predictions on test set

arguments = [
    '--test_path', datapath,
    '--preds_path', prediction_path,
    '--checkpoint_dir', ensemble_savedir,
    '--gpu', '0',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling'
]

args = chemprop.args.PredictArgs().parse_args(arguments)
preds = chemprop.train.make_predictions(args=args)

# %% Evaluates predictions by comparing to true values

