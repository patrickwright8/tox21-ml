# %% Evaluates PRC-AUC of predictions
import os
import pandas as pd
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

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
# %%
