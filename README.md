# Dyna_vivo_seq pipline

#### This pipline description mainly revolves around the generation of expression matrices with metabolically labeled transcripts using Dyna-vivo-seq, the generation of new and old transcript count matrices without labels, and the calculation of TC base mutation rates .

- ###### Main steps

[Step1](): Utilize the Drop-seq pipeline (James Nemesh, McCarroll Lab, version 1.12; Macosko et al., 2015) and the STAR v2.7.3a workflow for reads trimming and  alignment. In addition,both exonic and intronic regions that mapped to annotated gene loci were retained for downstream analysis.This step Modified from well-temp-seq.

```
python ~/Matrix_generation/drop-seq/Drop-seq_step1_GRCm38.py -F1
../raw/sample_name_S3_L006_R1_001.fastq.gz -F2 ../raw/sample_name_S3_L006_R2_001.fastq.gz -SM sample_name
```

The following steps refer to the scNT pipeline

[Step2]():Identification of T-to-C substitutions in control and experimental samples (without IAA treatment, as background)

```
sh ~/Matrix_generation/Step2_extract_alignment_info_GRCm38.sh . sample_name
```

[Step3]():exclude the genomic sites with background T-to-C substitutions .

```
###For samples without IAA treatment, a background deduction cannot be made using the following command
perl ~/Matrix_generation/drop-seq/scripts/TagIntronicRead_V5.pl -read
sample1_name_both_strand_all_TC.tsv_q27.tsv -bam ../sample1_name_star_gene_exon_tagged_TagIntronic_clean.bam
###The sample with IAA treatment needs to subtract background mutations through a control group. Use the following command
sh ~/Matrix_generation/Step3_substract_background_locus.sh sample_name2 sample_name1 ../ . ../../sample_name1/Step2/
```

[Step4]():Labeled and unlabeled transcripts gene expression matrix output

```
sh ~/Matrix_generation/Step4_genetare_TC_matrix.sh . sample_name cell_number
```

- ###### Unlabeled expression matrix extraction

```
python ~/Matrix_generation/Step2_DGE.py -b
*_star_gene_exon_tagged_TagIntronic_clean.bam -m F -n cell_number -sm sample_name
```

- ###### TC mutation proportion statistics

```
mkdir {cluster,label_ratio,library,mutation}
cd mutation
Rscript
~/plot/mutation_rate_cmd_v4.R  <input>  <file1>  <file2>  <file..> 
file123 : ../sample_name3/Step2/*_both_strand_all.tsv_q27_gene_anno_stat.txt
```

- Generate rds without TC labeling
- Total RDS with noType

```
##Barcode was extracted from rds files with mutation tags
cd ../
mkdir data_noType
cd data_noType/
ls ../sample_name/Step2/*rds | while read id
do Rscript ~/plot/extra_cb_from_TC_type_rds.R
$id 7000
done

ls ../sample_name/Step2/*gene_cell_UMI_read.txt | while read id
do sample=$(basename $id
'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_U
MI_read.txt')
cb_file=${sample}_TC_cb_7000.txt
Rscript
~/plot/Generate_T_C_matrix_noType_ne
w_old_total_4cmd_v2.R $id $cb_file cellnum $sample total
done
```

- New old RDS-NoType

```
##Generate intermediate file *corrected_gene_cell_UMI_read.txt

python
~/plot/select_new_old_RNA_fromCountF
ile.py -p1 ../sample_name/Step2/ -p2 . -s sample_name


###Generate rds files of old and new RNA
Rscript Generate_T_C_matrix_NoType_new_old_all_4cmd_v2.R <input>
<total_RNA_num_barcode_file> <num_core_barcode> <sample_name>
<data_type>

###### new
ls sample_name*read_new.txt | while read id
do sample=$(basename $id
'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_U
MI_read_new.txt')
cb_file=${sample}_TC_cb_7000.txt
Rscript
~/plot/Generate_T_C_matrix_noType_ne
w_old_total_4cmd_v2.R $id $cb_file cellnum $sample new
done
###### old
ls sample_name*read_old.txt | while read id
do sample=$(basename $id
'_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_U
MI_read_old.txt')
cb_file=${sample}_TC_cb_7000.txt
Rscript
~/plot/Generate_T_C_matrix_noType_ne
w_old_total_4cmd_v2.R $id $cb_file cellnum $sample old
done
```
Data:Raw data files are available at NCBI Gene Expression Omnibus (GEO):
