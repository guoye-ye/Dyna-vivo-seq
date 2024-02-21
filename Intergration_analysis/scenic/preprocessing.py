#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:27:52 2021

@author: quinn
"""

import numpy as np
import scanpy as sc
import sys

adata = sc.read_csv(sys.argv[1], first_column_names=True)
adata.var_names_make_unique()

mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
sc.pp.filter_genes(adata, min_cells=20)

import loompy as lp
row_attrs = {"Gene": np.array(adata.var_names),}
col_attrs = {"CellID": np.array(adata.obs_names),"nGene": np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(), "nUMI": np.array(np.sum(adata.X.transpose(),axis=0)).flatten(),}
lp.create(sys.argv[2],adata.X.transpose(),row_attrs, col_attrs)

