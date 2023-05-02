# tox21-ml

As many as 30% of candidate pharmaceuticals in clinical trials cannot be commercialized due to unexpected toxicity. Here, I train a directed message passing neural net [(Chemprop, Yang et al)] [1] on Tox21, a publicly available toxicity dataset, to predict toxicity of small molecules on a variety of different assays. This is a work in progress, so more will come soon.
<br><br>
## _hyperopt.py: 
Performs hyperparameter optimization on the training set
## clustering.py:
Uses UMAP to visualize training, test, and evaluation sets in chemical space. 