#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 20:35:51 2021

@author: quinn
"""

import loompy as lp
import pandas as pd
import sys

input = sys.argv[1]+'_SCENIC.loom'
lf = lp.connect(input, mode='r+', validate=False)
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID) 
lf.close()
output = sys.argv[1] + '_SCENIC.csv'
auc_mtx.to_csv(output, sep=',', header=True, index=True)

