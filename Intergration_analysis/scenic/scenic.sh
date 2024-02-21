expression_mtx=$1 #loom file converted from rds
sample_name=$2
s1_output=${expression_mtx}_step1.tsv #a tsv file

tf_list='/home/songjia/newdisk/quinn/reference/cisTarget_database/hs_hgnc_tfs.txt'
database='/home/songjia/newdisk/quinn/reference/cisTarget_database/hg19/hg19-500bp-upstream-10species.mc9nr.feather'
annotatations_file='/home/songjia/newdisk/quinn/reference/cisTarget_database/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
s2_output=${expression_mtx}_step2.csv #a csv file

s3_output=${expression_mtx}_SCENIC.loom

pyscenic grn --num_workers 20 --output ${s1_output} --method grnboost2 ${expression_mtx} ${tf_list}
pyscenic ctx ${s1_output}  ${database} --annotations_fname ${annotatations_file} --expression_mtx_fname ${expression_mtx} --mode "dask_multiprocessing" --output ${s2_output} --num_workers 20 --mask_dropouts 

pyscenic aucell ${expression_mtx}  ${s2_output} --output ${s3_output} --num_workers 20

