sample_name=$1

csvname_new=${sample_name}_new_raw.csv
loomname_new=${sample_name}_new_filter.loom
csvname_old=${sample_name}_old_raw.csv
loomname_old=${sample_name}_old_filter.loom

python /home/songjia/newdisk/quinn/scripts/scenic/preprocessing.py ${csvname_new} ${loomname_new}
python /home/songjia/newdisk/quinn/scripts/scenic/preprocessing.py ${csvname_old} ${loomname_old}


