# %% Imports, I/O

### Note: This program must be ran using the command line - I had trouble getting it to run using Spyder.

import chemprop
from chemprop.train.evaluate import evaluate_predictions
import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = os.getcwd()
datapath = os.path.join(root, 'data\\test_data\\combined_test_set.csv')
ensemble_savedir = os.path.join(root, 'ensemble_ckpts')
prediction_path = os.path.join(root, 'ensemble_preds\\predictions.csv')

# %% Run ensemble predictions on test set

arguments = [
    '--test_path', datapath,
    '--preds_path', prediction_path,
    '--checkpoint_dir', ensemble_savedir,
    '--num_workers', '0',
    '--gpu', '0',
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling',
    '--smiles_columns', 'smiles'
]

args = chemprop.args.PredictArgs().parse_args(arguments)
preds = chemprop.train.make_predictions(args=args)

# %% Evaluates predictions by comparing to true values

root = os.getcwd()
prediction_path = os.path.join(root, 'ensemble_preds\\predictions.csv')
formatted_datapath = os.path.join(root, 'data\\test_data\\formatted_combined_test_set.csv')

predictions = pd.read_csv(prediction_path)
targets = pd.read_csv(formatted_datapath)

predictions_as_list = predictions.iloc[:,1:].values.tolist()
targets_as_list = targets.iloc[:,1:].values.tolist()

# Data to plot precision - recall curve
precision, recall, thresholds = precision_recall_curve(targets_as_list, predictions_as_list)
# Use AUC function to calculate the area under the curve of precision recall curve
auc_precision_recall = auc(recall, precision)
print(auc_precision_recall)

'''
evaluate_predictions(
    preds=predictions_as_list,
    targets=targets_as_list,
    num_tasks=12,
    dataset_type='classification',
    metrics=['roc-auc']
)
'''
# %%
