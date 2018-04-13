# NGS_BRCA_Pipelines
Pipelines are developed to analyze NGS data for breast and lung cancers.

Steps to run the pipeline scripts:

1 trim the raw NGS read pairs;

2 align the trimmed fastq files to the reference DNA sequences;

3 perform the post-alignment filtration to discard low-quality alignment results and identify gene-specific primers within R1;

4 cluster the uniqe molecular barcodes within R2 and consolidate the read pairs having the same barcode sequences.

5 reformat the alignment sam file, converting the position-related columns to chromosome-wise locations.

6 call variants by employing smCounter.

7 annotate by clinvar and cosmic databases.
