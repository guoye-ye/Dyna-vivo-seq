sample_name=$1
filter_loom_new=${sample_name}_new_filter.loom
filter_loom_old=${sample_name}_old_filter.loom

sh /home/songjia/newdisk/quinn/scripts/scenic/preprocess.sh ${sample_name}

sh /home/songjia/newdisk/quinn/scripts/scenic/scenic.sh ${filter_loom_new}
sh /home/songjia/newdisk/quinn/scripts/scenic/scenic.sh ${filter_loom_old} 

python /home/songjia/newdisk/quinn/scripts/scenic/loomtocsv.py ${filter_loom_new}
python /home/songjia/newdisk/quinn/scripts/scenic/loomtocsv.py ${filter_loom_old}

