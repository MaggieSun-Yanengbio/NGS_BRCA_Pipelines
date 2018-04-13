#!/bin/sh
trim="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/trim_reads_qiagen_20180305.py"
cluster="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/cluster_barcode_gatk_20180328.py"
align="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/align_reads_gatk_20180329.py"
filter="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/post_alignment_filter_gatk_20180408.py"
variant_call="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/variant_call_gatk_20180328.py"
annotate="/home/administrator/pipeline/NGS_BRCA_Pipelines-master/annotation_gatk_20180330.py"

source="/home/administrator/data/gatk_0202/"
sample_names="T-1W_S5 T-1Z_S15 T-2W_S6 T-2Z_S16 T-3N1_S11 T-3N2_S12 T-3Q1_S13 T-3Q2_S14 T-3T1_S9 T-3T2_S10 U00_S4 U1_S3 Undetermined_S0 Z1_S2 ZQ_S1"
tail_name="L001"
out_dir="undetermined/"
bwa_dir="bwa"
index_name="/home/administrator/source/refSeq2"
ref_fa_file="/home/administrator/source/target_breast.refSeq2.fa"
samtools_dir="samtools"
primers="/home/administrator/source/DHS-001Z_primers_target.csv"
min_consolidate_qual=10
min_consolidate_freq=0.8
smcounter="/home/administrator/smCounter/smCounter-master/smCounter.py"
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
cosmic="/home/administrator/database/COSMIC_variant.csv"
clinvar="/home/administrator/database/breast_cancer_variant_clinvar.csv"
ref_ens="/home/administrator/source/geneid_cancer_qiagen.csv"
gatk="/home/administrator/gatk/gatk-package-4.0.1.2-local.jar"
knowkites1="/home/administrator/knowkites/1000G_phase1.indels.hg19.vcf"
knowkites2="/home/administrator/knowkites/1000G_phase1.snps.high_confidence.hg19.vcf"
knowkites3="/home/administrator/knowkites/Mills_and_1000G_gold_standard.indels.hg19.vcf"

for sample_name in $sample_names
do
  command_trim="$trim $source $sample_name $tail_name $out_dir"
  command_cluster="$cluster $source $sample_name $min_consolidate_qual $min_consolidate_freq"
  command_align="$align $source $sample_name $bwa_dir $index_name $ref_fa_file"
  command_post_alignment="$filter $source $sample_name $samtools_dir $primers"
  command_variant_call="$variant_call $source $sample_name $gatk $samtools_dir $refGenome $knowkites1 $knowkites2 $knowkites3"
  command_annotate="$annotate $source $sample_name $cosmic $clinvar $ref_ens"

  echo "0 Trim"
  python3 $command_trim
  echo "1 Clustering based on UMIs"
  python3 $command_cluster
  echo "2 Alignment"
  python3 $command_align
  echo "3 Post alignment analysis"
  python3 $command_post_alignment
  echo "4 Variant call"
  python3 $command_variant_call
  echo "5 Annotate"
  python3 $command_annotate
  echo "----------------------------------------------------------"
  echo "sample $sample_name is OK."
  echo "----------------------------------------------------------"
done