# NGS_BRCA_Pipelines
Pipelines are developed to analyze NGS data for breast and lung cancers.

Steps to run the pipeline scripts:

1 trim the raw NGS read pairs;

2 align the trimmed fastq files to the reference DNA sequences;

3 perform the post-alignment filtration to discard low-quality alignment results and identify gene-specific primers within R1;

4 cluster the uniqe molecular barcodes within R2 and consolidate the read pairs having the same barcode sequences;

5 reformat the alignment sam file, converting the position-related columns to chromosome-wise locations;

6 call variants by employing smCounter;

7 annotate by clinvar and cosmic databases.

------------------------------------------------------------------------------------------------------------

Run smCounter pipeline by run_pipeline.sh:

trim_reads.py: input: R1.fastq.gz, R2.fastq.gz
               output: R1_undetermined.fastq, R2_undetermined.fastq, basic_stats.txt

align_reads.py: input: R1_undetermined.fastq, R2_undetermined.fastq, refSeq.fa, refSeq index file(.amb .ann .bwt .pac .sa)
                output: aligned.sam, align_stats.txt
                
post_alignment_filter.py: input: aligned.sam, primers_file, samtools
                          output: filtered.sam, align_stats.txt

cluster_barcode.py: input: R1_undetermined.fastq, R2_undetermined.fastq
                    output: R1_consolidated.fastq, R2_consolidated.fastq, log, 
                    R1_undetermined.fastq_umi.fastq_ids.txt, R2_undetermined.fastq_umi.fastq_ids.txt

reformat_samfile.py: input: R1_consolidated.fastq, R2_consolidated.fastq, R1_undetermined.fastq_umi.fastq_ids.txt, filtered.sam
                     output: vcready.sam
                     
varaint_call.py: input: vcready.sam, samtools, bedtools, smcounter
                        "vcready.sam --- vcready.bam --- vcready_sorted.bam --- vcready_sorted.bam.bai" by samtools
                 output: smCounter.cut.vcf, log
                 
annotate.py: input: smCounter.cut.vcf, databases:cosmic clinvar g1000, reference Ensembl list:"genename,ENSG,ENST"
             output: annotated.csv, annotate_stats.csv
