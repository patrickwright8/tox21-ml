import pandas as pd
import numpy as np
import os

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
IPythonConsole.molSize = (350, 350)   # Change image size
IPythonConsole.ipython_useSVG = True  # Change output to SVG

root = os.getcwd() 
datapath = os.path.join(root, 'projects\\tox21-ml\\data\\test_data\\tox21_10k_challenge_test.sdf')
