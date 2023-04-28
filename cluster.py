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

# %% Combines train and eval sets into one df

train = pd.read_csv(train_path)
eval = pd.read_csv(eval_path)

train = train.iloc[:,3:]
eval = eval.iloc[:,3:]

train.insert(0, 'Type', 'Train')
eval.insert(0, 'Type', 'Eval')

df = pd.concat([train, eval])
df_vecs = df.iloc[:,1:]

# %% Trains UMAP on molecular vectors, then plots results

reducer = umap.UMAP(n_neighbors=15, min_dist=.05)

# Output of mol2vec is already "scaled", so no further input scaling is needed
embedding = reducer.fit_transform(df_vecs)

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=[sns.color_palette()[x] for x in df.Type.map({"Train":0, "Eval":1})])
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP projection of Mol2Vec Data', fontsize=24)
# %%
