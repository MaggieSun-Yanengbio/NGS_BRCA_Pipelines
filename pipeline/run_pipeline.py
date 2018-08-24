import os
import numpy as np
from collections import defaultdict

from qc import *
from trim import *
from alignment import *
from post_alignment_filter import *
from umi_cluster import *
from reformat_samfile import *
from variant_call import *
from annotation import *


# common parameters
source = '/home/administrator/pipeline_test/'
sample = 'N223-M2_S10'
tailname = 'L001'
stats_out = source + sample + 'statistics.csv'
read1 = source + sample + '_' + tailname + '_R1_001.fastq.gz'
read2 = source + sample + '_' + tailname + '_R2_001.fastq.gz'
fastqc_path = '/home/administrator/fastqc/FastQC/fastqc'
bwa_path = '/usr/bin/bwa'
samtools_path = 'samtools'
umitools_path = '/home/administrator/umitools/umi_tools/umi_tools.py'
bedtools_path = '/usr/bin/bedtools'
ref_genome = '/home/administrator/source/ucsc.hg19.fasta'
ref_fasta = '/home/administrator/source/target_breast.refSeq.fa'
ref_bed = '/home/administrator/source/target_breast.refSeq.bed'
index_name = '/home/administrator/source/refseq'
primers = '/home/administrator/source/DHS-001Z_primers_target.csv'
smcounter_path = '/home/administrator/smCounter/smCounter-master/smcounter_v1.py'
bed_tandem_repeats = '/home/administrator/smCounter/smCounter-master/simpleRepeat.bed'
bed_repeat_masker_subset = '/home/administrator/smCounter/smCounter-master/SR_LC_SL.nochr.bed'
geneinfo = '/home/administrator/source/geneid_breast_cancer.csv'
cosmic = '/home/administrator/database/cosmic_breast_93genes_20180820.csv'
clinvar = '/home/administrator/database/clinvar_breast_93genes_20180729.csv'
g1000 = '/home/administrator/database/1000genomes_breast_93genes.csv'


def main():
    # raw data qc parameters
    qc_dir = source + 'rawdata_qc'
        # raw data qc process
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    qc_result1, qc_result2 = qc_raw_reads(fastqc_path, qc_dir, sample, read1, read2,)
    min_length1, max_length1, gc1, q20_r1, q30_r1 = qc_result1
    min_length2, max_length2, gc2, q20_r2, q30_r2 = qc_result2
    raw_min_length = min(int(min_length1), int(min_length2))
    raw_max_length = max(int(max_length1), int(max_length2))
    raw_gc = round(np.mean([float(gc1), float(gc2)]), 3)
    raw_q20 = round(np.mean([float(q20_r1[:-1]), float(q30_r1[:-1])]), 3)
    raw_q30 = round(np.mean([float(q20_r2[:-1]), float(q30_r2[:-1])]), 3)


    # trim parameters
    out_dir = 'undetermined/'
    min_read_len = 40
    log_dir = source + out_dir + 'log/'
    trimmed1 = source + out_dir + sample + '_R1_undetermined.fastq'
    trimmed2 = source + out_dir + sample + '_R2_undetermined.fastq'
    trim_stats_file = source + out_dir + sample + '_basic_stats.txt'
    common_seq1 = 'CAAAACGCAATACTGTACATT'
    common_seq2 = 'ATTGGAGTCCT'
    # trim process
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    logger_trim_process, logger_trim_errors = store_logs(log_dir)
    num_total_reads, num_short_reads, num_unproper_reads, num_without_comseq_reads, num_clean_reads, percentage_clean = \
        trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len, common_seq1, common_seq2,
                        trim_stats_file, logger_trim_process, logger_trim_errors)


    # clean data qc parameters
    qc_dir = source + 'cleandata_qc'
    # clean data qc process
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    qc_result1, qc_result2 = qc_raw_reads(fastqc_path, qc_dir, sample, trimmed1, trimmed2)
    min_length1, max_length1, gc1, q20_r1, q30_r1 = qc_result1
    min_length2, max_length2, gc2, q20_r2, q30_r2 = qc_result2
    clean_min_length = min(int(min_length1), int(min_length2))
    clean_max_length = max(int(max_length1), int(max_length2))
    clean_gc = round(np.mean([float(gc1), float(gc2)]), 3)
    clean_q20 = round(np.mean([float(q20_r1[:-1]), float(q30_r1[:-1])]), 3)
    clean_q30 = round(np.mean([float(q20_r2[:-1]), float(q30_r2[:-1])]), 3)


    # parameter alignment
    out_file = source + 'aligned/' + sample + '_aligned.sam'
    log_dir = source + 'aligned/log/'
    num_threads = 4
    # alignment process
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    length_target_region = 0
    for row in open(ref_fasta, 'r'):
        if row.startswith('>'):
            start, end = row.strip().split('_')[1:]
            length_target_region += int(end) - int(start)
    logger_bwa_process, logger_bwa_errors = store_logs(log_dir)
    align_reads_bwa(bwa_path, ref_fasta, index_name, trimmed1, trimmed2,
                    out_file, num_threads, logger_bwa_process, logger_bwa_errors)


    # filter parameters
    log_dir = source + 'aligned/log/'
    alignment_sam = source + 'aligned/' + sample + '_aligned.sam'
    min_mapq = 17
    max_soft_clip = 10
    out_file1 = source + 'aligned/' + sample + '_tmp.sam'
    filter_stats_file = source + 'aligned/' + sample + '_align_stats.txt'
    primer_stats_file = source + 'aligned/' + sample + '_primer_stats.csv'
    max_dist = 2
    out_file = source + 'aligned/' + sample + '_filtered.sam'
    # filter process
    logger_filter_process, logger_filter_errors = store_logs(log_dir)
    (num_sec, num_unmap, min_mapq, num_low_mapq, final_count, percentage_filter_reads) = \
        filter_alignment_samtools(samtools_path, alignment_sam, min_mapq,
                                max_soft_clip, out_file1, filter_stats_file,
                                logger_filter_process, logger_filter_errors)
    (num_unproper_pairs, num_primers, num_off_target, ratio_off, num_on_target, ratio_on, num_target_primers, mean_rdp,
     num_primers_5mean_rdp, num_primers_25mean_rdp, num_primers_50mean_rdp, num_primers_75mean_rdp, num_primers_mean_rdp,
     mean_mdp, num_primers_5mean_mdp, num_primers_25mean_mdp, num_primers_50mean_mdp, num_primers_75mean_mdp, num_primers_mean_mdp
     )= identify_gs_primers(samtools_path, out_file1, primers, max_dist, out_file,
                            filter_stats_file, primer_stats_file, logger_filter_process,
                            logger_filter_errors)

    # umi cluster parameters
    filein_sam = source + 'aligned/' + sample + '_filtered.sam'
    filein_bam = source + 'aligned/' + sample + '_filtered.bam'
    sorted_bam = source + 'aligned/' + sample + '_filtered_sorted.bam'
    clustered_bam = source + 'aligned/' + sample + '_clustered.bam'
    clustered_sam = source + 'aligned/' + sample + '_clustered.sam'
    stats = source + 'aligned/' + sample + '_deduplicated'
    coverage_file = source + 'aligned/' + sample + '_coverage.tsv'
    # umi cluster process
    (total_num_bases, mean_bd, min_bd, max_bd, num_50x, num_100x, num_200x, num_500x,
     num_5mean, num_25mean, num_50mean, num_75mean, num_mean) \
        = stats_target_bases(samtools_path, sorted_bam, ref_bed, coverage_file)
    umi_depth_to_cnt = defaultdict(int)
    depth_distribution = [0]
    umitools_return = umitools(samtools_path, umitools_path, filein_sam, filein_bam, sorted_bam, clustered_bam, clustered_sam, stats)
    for row in umitools_return.decode('utf-8').split('\n'):
        if not row.startswith('#'):
            if '#umis' in row:
                num_mts = row.split(' ')[-1].strip()
            if 'Input Reads' in row:
                num_input = int(row.split(' ')[-1].strip())
            if 'reads out' in row:
                num_output = int(row.split(' ')[-1].strip())
    num_consilidated = num_input - num_output

    per_umi_stats = open(stats + '_per_umi.tsv','r')
    per_umi_stats.readline()
    for row in per_umi_stats.readlines():
        depth = row.strip().split('\t')[3]
        umi_depth_to_cnt[int(depth)] += 1
    for dp in sorted(umi_depth_to_cnt.keys()):
        depth_distribution += [dp] * umi_depth_to_cnt[dp]
    mt_depth_25 = int(np.percentile(depth_distribution, 25))
    mt_depth_median = int(np.percentile(depth_distribution, 50))
    mt_depth_75 = int(np.percentile(depth_distribution, 75))
    mt_depth_max = int(np.percentile(depth_distribution, 100))

    # reformat samfile parameters
    alignment_sam = source + 'aligned/' + sample + '_clustered.sam'
    vcready_sam = source + 'aligned/' + sample + '_vcready.sam'
    # reformat samfile process
    reformat_sam(alignment_sam, vcready_sam)
    vcready_bam = sam_to_bam(vcready_sam, samtools_path)
    vcready_sorted_bam = sort_index(vcready_bam, samtools_path)

    # variant calling parameters
    normal = '/home/administrator/pipeline_test/aligned/N223-N_S9_vcready_sorted.bam'
    outprefix = source + sample + '_variant'
    logfile = source + sample + '_logfile'
    n_cpu = '16'
    min_bq = '20'
    min_mq = '20'
    rpb = '8.6'
    mt_drop = '0'
    hp_len = '8'
    mismatch = '6'
    threshold = '0'
    min_frequency = '0.01'
    min_active_score = '0.01'
    tlod_threshold = '6.3'
    nlod_threshold = '3.5'
    # variant calling process
    call(smcounter_path, outprefix, vcready_sorted_bam, normal, ref_bed, ref_fasta, rpb, n_cpu, min_bq, min_mq, hp_len,
         mismatch, mt_drop, threshold, min_frequency, min_active_score, tlod_threshold, nlod_threshold, ref_genome,
         bed_tandem_repeats, bed_repeat_masker_subset, bedtools_path, logfile)

    # annotation parameters
    variant_vcf = source + sample + '_variant.smCounter.somatic.cut.vcf'
    annotated_csv = source + sample + '_annotated.csv'
    annotation_stats_file = source + sample + '_annotate_stats.txt'
    # annotation process
    dict_cos, dict_clin, dict_g1000 = read_database(cosmic,clinvar,g1000)
    annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,annotation_stats_file)
    fill_table(annotated_csv, geneinfo)


    # statistics parameters
    pipeline_stats_file = source + sample + '_statistics_in_report.csv'
    # statistics process
    with open(pipeline_stats_file, 'w') as f:
        f.write('library name,' + sample + '\n')
        f.write('Number of raw read pairs,' + str(num_total_reads) + '\n')
        f.write('Minimun read length (raw),' + str(raw_min_length) + '\n')
        f.write('Maximun read length (raw),' + str(raw_max_length) + '\n')
        f.write('GC%(raw),' + str(raw_gc) + '\n')
        f.write('Q20 in raw data(%),' + str(raw_q20) + '\n')
        f.write('Q30 in raw data(%),' + str(raw_q30) + '\n')
        f.write('Number of short read pairs (either read_length < 40bp),' + str(num_short_reads) + '\n')
        f.write('Number of unproper read pairs (containing incorrect headers),' + str(num_unproper_reads) + '\n')
        f.write('Number of read pairs without correct common sequences/Mts,' + str(num_without_comseq_reads) + '\n')
        f.write('Number of clean read pairs,' + str(num_clean_reads) + '\n')
        f.write('Percentage of clean read pairs(%),' + str(percentage_clean) + '\n')
        f.write('Minimun read length (clean),' + str(clean_min_length) + '\n')
        f.write('Maxnimun read length (clean),' + str(clean_max_length) + '\n')
        f.write('GC%(clean),' + str(clean_gc) + '\n')
        f.write('Q20 in clean data(%),' + str(clean_q20) + '\n')
        f.write('Q30 in clean data(%),' + str(clean_q30) + '\n')
        f.write('Number of secondary/supplimentary alignments,' + str(num_sec) + '\n')
        f.write('Number of unmapped read pairs,' + str(num_unmap) + '\n')
        f.write('Threshold MAPQ for filtration,' + str(min_mapq) + '\n')
        f.write('Number of low MAPQ alignments(MAPQ<threshold),' + str(num_low_mapq) + '\n')
        f.write('Number of alignments passed SAMTools filtration,' + str(final_count) + '\n')
        f.write('Percentage of alignments passed SAMTools filtration(%),' + str(percentage_filter_reads) + '\n')
        f.write('Number of unproperly-paired alignments,' + str(num_unproper_pairs) + '\n')
        f.write('Number of panel primers,' + str(num_primers) + '\n')
        f.write('Number of off-target alignments,' + str(num_off_target) + '\n')
        f.write('Percentage of off-target alignments(%),' + str(ratio_off) + '\n')
        f.write('Number of on-target alignments,' + str(num_on_target) + '\n')
        f.write('Percentage of on-target alignments(%),' + str(ratio_on) + '\n')
        f.write('Number of MTs,' + num_mts + '\n')
        f.write('Number of consolidated read pairs,' + str(num_consilidated) + '\n')
        f.write('Mean read depth of MT,' + str(round(num_on_target/int(num_mts), 3)) + '\n')
        f.write('25% of read depth per MT,' + str(mt_depth_25) + '\n')
        f.write('Median read depth per MT,' + str(mt_depth_median) + '\n')
        f.write('75% of read depth per MT,' + str(mt_depth_75) + '\n')
        f.write('Maximun read depth per MT,' + str(mt_depth_max) + '\n')
        f.write('Number of targeted primers,' + str(num_target_primers) + '\n')
        f.write('Mean MT depth per primer (mean MDP),' + str(mean_mdp) + '\n')
        f.write('% of primers >= 5% of mean MDP,' + str(round(100 * num_primers_5mean_mdp / num_primers, 3)) + '\n')
        f.write('% of primers >= 25% of mean MDP,' + str(round(100 * num_primers_25mean_mdp / num_primers, 3)) + '\n')
        f.write('% of primers >= 50% of mean MDP,' + str(round(100 * num_primers_50mean_mdp / num_primers, 3)) + '\n')
        f.write('% of primers >= 75% of mean MDP,' + str(round(100 * num_primers_75mean_mdp / num_primers, 3)) + '\n')
        f.write('% of primers >= mean MDP,' + str(round(100 * num_primers_mean_mdp / num_primers, 3)) + '\n')
        f.write('Mean read depth per primer (mean RDP),' + str(mean_rdp) + '\n')
        f.write('% of primers >= 5% of mean RDP,' + str(round(100 * num_primers_5mean_rdp/num_primers, 3)) + '\n')
        f.write('% of primers >= 25% of mean RDP,' + str(round(100 * num_primers_25mean_rdp/num_primers, 3)) + '\n')
        f.write('% of primers >= 50% of mean RDP,' + str(round(100 * num_primers_50mean_rdp/num_primers, 3)) + '\n')
        f.write('% of primers >= 75% of mean RDP,' + str(round(100 * num_primers_75mean_rdp/num_primers, 3)) + '\n')
        f.write('% of primers >= mean RDP,' + str(round(100 * num_primers_mean_rdp/num_primers, 3)) + '\n')
        f.write('Length of target region,' + str(length_target_region) + '\n')
        f.write('Total number of sequencing bases (covering target region),' + str(total_num_bases) + '\n')
        f.write('Mean sequencing depth per base (mean BD),' + str(mean_bd) + '\n')
        f.write('Minimun BD,' + str(min_bd) + '\n')
        f.write('Maxnimun BD,' + str(max_bd) + '\n')
        f.write('% of target bases with BD >= 50x,' + str(num_50x) + '\n')
        f.write('% of target bases with BD >= 100x,' + str(num_100x) + '\n')
        f.write('% of target bases with BD >= 200x,' + str(num_200x) + '\n')
        f.write('% of target bases with BD >= 500x,' + str(num_500x) + '\n')
        f.write('% of target bases >= 5% of mean BD,' + str(num_5mean) + '\n')
        f.write('% of target bases >= 25% of mean BD,' + str(num_25mean) + '\n')
        f.write('% of target bases >= 50% of mean BD,' + str(num_50mean) + '\n')
        f.write('% of target bases >= 75% of mean BD,' + str(num_75mean) + '\n')
        f.write('% of target bases >= mean BD,' + str(num_mean) + '\n')


if __name__ == '__main__':
    main()
