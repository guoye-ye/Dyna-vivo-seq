# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 16:24:03 2020

@author: cyclopenta

"""
#### v2和v1的区别就是v2可以通过指定用于分析标记比例的细胞数目

from collections import defaultdict, Counter
import re, pandas
import sys
umi_gene_count = defaultdict(Counter)
umi_gene_snp = defaultdict(dict)
umi_gene_assign = {}
barcode_gene_snp = defaultdict(Counter)
refined_mutation_dict = defaultdict(Counter)
snp_bool = 0
tmp_c = 0
file = open(sys.argv[1])#*tidy.txt
file2 = open(sys.argv[2])#*cb.txt
num = int(sys.argv[3])#cell_barcode

for line in file:
    tmp_c += 1
    if tmp_c == 1:
        continue
    line = line.rstrip()
    info_list = line.split(' ')
    gene = info_list[0]
    snp = info_list[1]
    umi_ubc = info_list[2] + info_list[3]
    umi_counts = int(info_list[4])
    umi_gene_count[umi_ubc][gene] += umi_counts
    umi_gene_snp[umi_ubc][gene] = snp
umi_gene_count= dict(umi_gene_count)
umi_gene_snp = dict(umi_gene_snp)

for k, v in umi_gene_count.items():
    umi_gene = max(v,key = v.get)
    umi_gene_assign[k] = umi_gene

for k, v in umi_gene_assign.items():
    barcode = k[0:12]
    snp_type = umi_gene_snp.get(k, 'error').get(v, 'error')
    if snp_type == 'T':
        barcode_gene_snp[barcode][snp_type] += 1
    elif snp_type =='C':
        barcode_gene_snp[barcode][snp_type] += 1
    else:
        print('error: ', v, k)
barcode_gene_snp = dict(barcode_gene_snp)

'''
edit
'''


cell_barcode = file2.readlines()
i = 0
for i in range(0,num):
    cell_barcode[i] = cell_barcode[i].rstrip()

for k in barcode_gene_snp.keys():
    if k in cell_barcode:
        mutation = barcode_gene_snp[k]
        refined_mutation_dict[k] = mutation
Data = pandas.DataFrame(refined_mutation_dict).T.fillna(0)
Data.to_csv(sys.argv[4])#export txt name

'''for v in barcode_gene_snp.values():
    for v2 in v.values():
        if v2 >= 2:
            print(v)'''

