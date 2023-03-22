# %% Imports, I/O
import chemprop
import pandas as pd
import numpy as np 
import pickle as pkl
import os

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = 'C:/projects/tox21-ml' # Eventually change to os.getcwd()
datafolder = root + '/data/'
datapath = datafolder + 'combined_tox21.pkl'

with open(datapath, 'rb') as r:
    df = pkl.load(r)
    
    