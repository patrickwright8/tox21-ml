# %% Imports, I/O
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import umap

sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})

root = os.getcwd()
train_path = os.path.join(root, 'data\\train2vec.csv')
eval_path = os.path.join(root, 'data\\test_data\\eval2vec.csv')
test_path = os.path.join(root, 'data\\test_data\\test2vec.csv')

# %% Combines train, test, and eval sets into one df

train = pd.read_csv(train_path)
eval = pd.read_csv(eval_path)
test = pd.read_csv(test_path)

train = train.iloc[:,3:]
eval = eval.iloc[:,3:]
test = test.iloc[:,3:]

train.insert(0, 'Type', 'Train')
eval.insert(0, 'Type', 'Eval')
test.insert(0, 'Type', 'Test')

df = pd.concat([train, test, eval])
df_vecs = df.iloc[:,1:]

# %% Trains UMAP on molecular vectors, then plots results

reducer = umap.UMAP()

# Output of mol2vec is already "scaled", so no further input scaling is needed
embedding = reducer.fit_transform(df_vecs)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=[sns.color_palette()[x] for x in df.Type.map({"Train":0, "Test":1, "Eval":2})])
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP projection of Mol2Vec Data', fontsize=24)
plt.legend(labels=['Train', 'Test', 'Eval'])
# %%
