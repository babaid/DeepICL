#!/bin/bash

# 1. Create conda environment
mamba create -n DeepICL python=3.7.13 -y
mamba activate DeepICL

# 2. install packages dependent on mkl
mamba install scipy=1.7.3 numpy=1.21.1 pandas=1.3.4 scikit-learn=1.0.2 seaborn=0.11.0 -y

# 3. install pytorch and torch_geometric
mamba install pytorch_geometric==2.0.3 -y
mamba install pytorch_cluster -y

# 4. install others
mamba install -c rdkit rdkit=2022.09.1 -y
mamba install -c conda-forge biopython=1.77 openbabel=3.1.1 -y
mamba install -c conda-forge plip=2.2.2


