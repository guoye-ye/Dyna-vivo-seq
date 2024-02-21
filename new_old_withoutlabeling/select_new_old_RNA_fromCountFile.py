#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 18:21:25 2021

@author: april
"""

"""
用最基础的count文件提取新/旧RNA，下一步再生成表达矩阵 (下一步使用)
I1:Jurkat_4sU_IAA_paired_seq_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read.txt
O1:Jurkat_4sU_IAA_paired_seq_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read_old.txt
O2:Jurkat_4sU_IAA_paired_seq_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read_new.txt
"""
import sys
import os
import argparse
import re
from collections import defaultdict

my_parse = argparse.ArgumentParser(description='rate of three type marks for specific gene')
my_parse.add_argument('-p1', '--PATH1', type = str, required = True,
                      help = 'the directory for raw_count.txt ') #eg, path/to/you
my_parse.add_argument('-p2', '--PATH2', type = str, required = True,
                      help = 'the directory where the save file in ') #eg, path/to/you
my_parse.add_argument('-s', '--SAMPLE', type = str, nargs='+', required=True, 
                      help='input the sample name') #输入几个样品的名称

args = my_parse.parse_args()

os.chdir('%s' % args.PATH2)


def selecet_newRNA():
    sample_list = args.SAMPLE
    infile1_list = []
    outfile1_list = []
    newRNA = []
    file_count = 0
    for filename in sample_list:
        filename = filename.rstrip()
        infile1 = str('%s'%args.PATH1)+'/'+str(filename)+'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read.txt'
        outfile1 = str(filename)+'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read_new.txt'
        #infile1 = str('%s'%args.PATH1)+'/'+str(filename)+'.txt'
        #outfile1 = str(filename)+'_new.txt'
        infile1_list.append(infile1)
        outfile1_list.append(outfile1)#得到输出文件名的list
    for infile1 in infile1_list:
        newRNA = []
        file_count += 1
        infile1 = infile1.rstrip()
        with open(infile1, 'r') as if1:
            for line in if1:
                line = line.rstrip()
                pattern1 = re.compile('--C\s')
                flag1 = re.search(pattern1, line)
                if flag1:
                #if line.find(str('%s'%args.GENE)+'--')>=0:
                    newRNA.append(line)
            print(str(sample_list[file_count-1])+'_the_newRNA_count:'+str(len(newRNA)))
        outfileName = outfile1_list[file_count-1]
        with open(outfileName, 'w') as of1:
            for line in newRNA:
                line = line.rstrip()
                of1.write(str(line)+'\n')

def selecet_oldRNA():
    sample_list = args.SAMPLE
    infile1_list = []
    outfile1_list = []
    oldRNA = []
    file_count = 0
    for filename in sample_list:
        filename = filename.rstrip()
        infile1 = str('%s'%args.PATH1)+'/'+str(filename)+'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read.txt'
        outfile1 = str(filename)+'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read_old.txt'
        #infile1 = str('%s'%args.PATH1)+'/'+str(filename)+'.txt'
        #outfile1 = str(filename)+'_old.txt'
        infile1_list.append(infile1)
        outfile1_list.append(outfile1)#得到输出文件名的list
    for infile1 in infile1_list:
        oldRNA = []
        file_count += 1
        infile1 = infile1.rstrip()
        with open(infile1, 'r') as if1:
            for line in if1:
                line = line.rstrip()
                pattern1 = re.compile('--T\s')
                flag1 = re.search(pattern1, line)
                if flag1:
                #if line.find(str('%s'%args.GENE)+'--')>=0:
                    oldRNA.append(line)
            print(str(sample_list[file_count-1])+'_the_oldRNA_count:'+str(len(oldRNA)))
        outfileName = outfile1_list[file_count-1]
        with open(outfileName, 'w') as of1:
            for line in oldRNA:
                line = line.rstrip()
                of1.write(str(line)+'\n')


selecet_newRNA()
selecet_oldRNA()