# tox21-ml

As many as 30% of candidate pharmaceuticals in clinical trials cannot be commercialized due to unexpected toxicity. Here, I train a directed message passing neural net [(Chemprop, Yang et al)] on Tox21, a publicly available toxicity dataset, to predict toxicity of small molecules on a variety of different assays. This is a work in progress, so more will come soon.
<br>
## _hyperopt.py: 
Performs hyperparameter optimization on the training set.
## clustering.py:
Applies UMAP to visualize training, test, and evaluation sets in chemical space. The input CSV files for this consist of the training, test, and evaluation datasets converted to vectors by Mol2vec. 

[(Chemprop, Yang et al)]: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00237