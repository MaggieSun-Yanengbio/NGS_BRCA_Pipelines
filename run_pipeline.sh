#!/bin/sh
trim="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/trim_reads_qiagen_20180305.py"
align="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/align_reads_qiagen_20180220.py"
filter="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/post_alignment_filter_qiagen_20180220.py"
cluster="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/cluster_barcode_qiagen_20180320.py"
reformat="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/reformat_samfile_qiagen_20180320.py"
variant_call="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/variant_call_qiagen_20180315.py"

source="/home/administrator/sample/"
sample_names="S0721_05B_CHG029767-YN-171205-N701-TAAGGCGA S0721_05B_CHG029767-YN-171205-N703-AGGCAGAA S0721_05B_CHG029767-YN-171205-N707-CTCTCTAC S0721_05B_CHG029767-YN-171205-N711-AAGAGGCA S0721_05B_CHG029767-YN-171205-N715-ATCTCAGG"
tail_name="L003"
out_dir="trimmed_"
bwa_dir="bwa"
index_name="/home/administrator/source/index"
ref_fa_file="/home/administrator/source/target_breast.refSeq.fa"
samtools_dir="samtools"
primers="/home/administrator/source/DHS-001Z_primers_target.csv"
min_consolidate_qual=10
min_consolidate_freq=0.8
smcounter="/home/administrator/smCounter/smCounter-master/smCounter.py"
header="/home/administrator/source/header.sam"
bedTarget="/home/administrator/source/DHS-001Z.covered-150bp.bed"
mtDepth="3612"
rpb="8.6"
nCPU="16"
minBQ="20"
minMQ="30"
hpLen="8"
mismatchThr="6"
mtDrop="1"
threshold="60"
refGenome="/home/administrator/source/hg19.fa"
bedTandemRepeats="/home/administrator/smCounter/smCounter-master/simpleRepeat.bed"
bedRepeatMaskerSubset="/home/administrator/smCounter/smCounter-master/SR_LC_SL.nochr.bed"
bedtoolsPath="/usr/bin/"

for sample_name in $sample_names
do
  command_trim="$trim $source $sample_name $tail_name $out_dir"
  command_align="$align $source $sample_name $bwa_dir $index_name $ref_fa_file"
  command_post_alignment="$filter $source $sample_name $samtools_dir $primers"
  command_cluster="$cluster $source $sample_name $min_consolidate_qual $min_consolidate_freq"
  command_reformat_sam="$reformat $source $sample_name"
  command_variant_call="$variant_call $source $sample_name $samtools_dir $smcounter $bedTarget $mtDepth $rpb $nCPU $minBQ $minMQ $hpLen $mismatchThr $mtDrop $threshold $refGenome $bedTandemRepeats $bedRepeatMaskerSubset $bedtoolsPath"

  echo "0 Trim"
  python3 $command_trim
  echo "1 Alignment"
  python3 $command_align
  echo "2 Post alignment analysis"
  python3 $command_post_alignment
  echo "3 Clustering based on UMIs"
  python3 $command_cluster
  echo "4 Reformat SAM files"
  python3 $command_reformat_sam
  echo "5 Variant call"
  python3 $command_variant_call
done
