#!/bin/sh
source_data="/home/yaneng/RSun/Data/NGS2018_02_02/"
sample_names="ZQ_S1 T-1W_S5 T-1Z_S15"
bwa_dir="/home/yaneng/RSun/Softwares/bwa/bwa"
index_name="/home/yaneng/RSun/Softwares/bwa/QIAGEN_DHS_001Z"
ref_fa_file="/home/yaneng/RSun/Data/qiagen-breast/target_breast.refSeq.fa"
samtools_dir="/home/yaneng/RSun/Softwares/samtools/samtools"
primers="/home/yaneng/RSun/Data/qiagen-breast/DHS-001Z_primers_target.csv"
min_consolidate_qual=10
min_consolidate_freq=0.8

for sample_name in $sample_names
do
  command_align="/home/yaneng/RSun/Projects/qiagen-seq/pipelines/version02/align_reads_qiagen_20180220.py $source_data $sample_name $bwa_dir $index_name $ref_fa_file"
  command_post_alignment="/home/yaneng/RSun/Projects/qiagen-seq/pipelines/version02/post_alignment_filter_qiagen_20180220.py $source_data $sample_name $samtools_dir $primers"
  command_cluster="/home/yaneng/RSun/Projects/qiagen-seq/pipelines/version02/cluster_barcode_qiagen_20180320.py $source_data $sample_name $min_consolidate_qual $min_consolidate_freq"
  command_reformat_sam="/home/yaneng/RSun/Projects/qiagen-seq/pipelines/version02/reformat_samfile_qiagen_20180320.py $source_data $sample_name"

  echo "1 Alignment"
  python3 $command_align
  echo "2 Post alignment analysis"
  python3 $command_post_alignment
  echo "3 Clustering based on UMIs"
  python3 $command_cluster
  echo "4 Reformat SAM files"
  python3 $command_reformat_sam
done 