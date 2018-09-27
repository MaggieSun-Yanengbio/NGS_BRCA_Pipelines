import os
import numpy as np
from collections import defaultdict
import yaml
import time

from qc import *
from trim import *
from alignment import *
from post_alignment_filter import *
from umi_cluster import *
from reformat_samfile import *
from variant_call import *
from annotation import *


def main():

    time_start = time.time()

    yamlfile = '/home/administrator/pipeline/NGS_BRCA_Pipelines-master/parameters.yml'
    y = yaml.load(open(yamlfile))

    print('#################################')
    print('1.raw data qc is starting')
    print('#################################')

    # raw data qc parameters
    qc_dir = y['source'] + 'rawdata_qc/'
    qc_log_dir = qc_dir + 'log/'
    read1 = y['source'] + y['sample'] + '_' + y['tailname'] + '_R1_001.fastq.gz'
    read2 = y['source'] + y['sample'] + '_' + y['tailname'] + '_R2_001.fastq.gz'
    # raw data qc process
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    if not os.path.exists(qc_log_dir):
        os.makedirs(qc_log_dir)

    logger_qc_process, logger_qc_errors = store_logs_qc(qc_log_dir)
    qc_result1, qc_result2 = qc_raw_reads(y['fastqc_path'], qc_dir, y['sample'], read1, read2, logger_qc_process, logger_qc_errors)
    min_length1, max_length1, gc1, q20_r1, q30_r1 = qc_result1
    min_length2, max_length2, gc2, q20_r2, q30_r2 = qc_result2
    raw_min_length = min(int(min_length1), int(min_length2))
    raw_max_length = max(int(max_length1), int(max_length2))
    raw_gc = round(np.mean([float(gc1), float(gc2)]), 3)
    raw_q20 = round(np.mean([float(q20_r1[:-1]), float(q30_r1[:-1])]), 3)
    raw_q30 = round(np.mean([float(q20_r2[:-1]), float(q30_r2[:-1])]), 3)
    print ('Min length of raw reads == {0}'.format(raw_min_length))
    print ('Max length of raw reads == {0}'.format(raw_max_length))
    print ('Percentage of GC of raw reads == {0}'.format(raw_gc))
    print ('Percentage of Q20 of raw reads == {0}'.format(raw_q20))
    print ('Percentage of Q30 of raw reads == {0}'.format(raw_q30))
    print ('Finish.')

    print('#################################')
    print('2.trim is starting')
    print('#################################')

    # trim parameters
    trim_dir = 'undetermined/'
    trim_log_dir = y['source'] + trim_dir + 'log/'
    trimmed1 = y['source'] + trim_dir + y['sample'] + '_R1_undetermined.fastq'
    trimmed2 = y['source'] + trim_dir + y['sample'] + '_R2_undetermined.fastq'
    trim_stats_file = y['source'] + trim_dir + y['sample'] + '_basic_stats.txt'
    common_seq1 = 'CAAAACGCAATACTGTACATT'
    common_seq2 = 'ATTGGAGTCCT'

    # trim process
    if not os.path.exists(trim_log_dir):
        os.makedirs(trim_log_dir)
    logger_trim_process, logger_trim_errors = store_logs_trim(trim_log_dir)
    num_total_reads, num_short_reads, num_unproper_reads, num_without_comseq_reads, num_clean_reads, percentage_clean = \
        trim_read_pairs(read1, read2, trimmed1, trimmed2, y['min_read_len'], common_seq1, common_seq2,
                        trim_stats_file, logger_trim_process, logger_trim_errors)
    print ('Finish.')

    print('#################################')
    print('3.clean data qc is starting')
    print('#################################')

    # clean data qc parameters
    qc_dir = y['source'] + 'cleandata_qc/'
    qc_log_dir = qc_dir + 'log/'

    # clean data qc process
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    if not os.path.exists(qc_log_dir):
        os.makedirs(qc_log_dir)
    logger_qc_process, logger_qc_errors = store_logs_qc(qc_log_dir)
    qc_result1, qc_result2 = qc_raw_reads(y['fastqc_path'], qc_dir, y['sample'], trimmed1, trimmed2, logger_qc_process, logger_qc_errors)
    min_length1, max_length1, gc1, q20_r1, q30_r1 = qc_result1
    min_length2, max_length2, gc2, q20_r2, q30_r2 = qc_result2
    clean_min_length = min(int(min_length1), int(min_length2))
    clean_max_length = max(int(max_length1), int(max_length2))
    clean_gc = round(np.mean([float(gc1), float(gc2)]), 3)
    clean_q20 = round(np.mean([float(q20_r1[:-1]), float(q30_r1[:-1])]), 3)
    clean_q30 = round(np.mean([float(q20_r2[:-1]), float(q30_r2[:-1])]), 3)
    print('Min length of clean reads == {0}'.format(clean_min_length))
    print('Max length of clean reads == {0}'.format(clean_max_length))
    print('Percentage of GC of clean reads == {0}'.format(clean_gc))
    print('Percentage of Q20 of clean reads == {0}'.format(clean_q20))
    print('Percentage of Q30 of clean reads == {0}'.format(clean_q30))
    print ('Finish.')

    print('#################################')
    print('4.alignment is starting')
    print('#################################')

    # parameter alignment
    alignment_sam = y['source'] + 'aligned/' + y['sample'] + '_aligned.sam'
    align_log_dir = y['source'] + 'aligned/log/'

    # alignment process
    if not os.path.exists(align_log_dir):
        os.makedirs(align_log_dir)
    length_target_region = 0
    for row in open(y['ref_fasta'], 'r'):
        if row.startswith('>'):
            start, end = row.strip().split('_')[1:]
            length_target_region += int(end) - int(start)
    logger_bwa_process, logger_bwa_errors = store_logs_align(align_log_dir)
    align_reads_bwa(y['bwa_path'], y['ref_fasta'], y['index_name'], trimmed1, trimmed2,
                    alignment_sam, str(y['num_threads']), logger_bwa_process, logger_bwa_errors)
    print ('Finish.')

    print('#################################')
    print('5.filter is starting')
    print('#################################')

    # filter parameters
    filter_log_dir = y['source'] + 'aligned/log/'
    tmp_sam = y['source'] + 'aligned/' + y['sample'] + '_tmp.sam'
    filter_stats_file = y['source'] + 'aligned/' + y['sample'] + '_align_stats.txt'
    primer_stats_file = y['source'] + 'aligned/' + y['sample'] + '_primer_stats.csv'
    filter_sam = y['source'] + 'aligned/' + y['sample'] + '_filtered.sam'

    # filter process
    logger_filter_process, logger_filter_errors = store_logs_filter(filter_log_dir)
    (num_sec, num_unmap, min_mapq, num_low_mapq, final_count, percentage_filter_reads) = \
        filter_alignment_samtools(y['samtools_path'], alignment_sam, y['min_mapq'],
                                y['max_soft_clip'], tmp_sam, filter_stats_file,
                                logger_filter_process, logger_filter_errors)
    (num_unproper_pairs, num_primers, num_off_target, ratio_off, num_on_target, ratio_on, num_target_primers, mean_rdp,
     num_primers_5mean_rdp, num_primers_25mean_rdp, num_primers_50mean_rdp, num_primers_75mean_rdp, num_primers_mean_rdp,
     mean_mdp, num_primers_5mean_mdp, num_primers_25mean_mdp, num_primers_50mean_mdp, num_primers_75mean_mdp, num_primers_mean_mdp
     )= identify_gs_primers(y['samtools_path'], tmp_sam, y['primers'], y['max_dist'], filter_sam,
                            filter_stats_file, primer_stats_file, logger_filter_process,
                            logger_filter_errors)
    print ('Finish.')

    print('#################################')
    print('6.umi cluster is starting')
    print('#################################')

    # umi cluster parameters
    umitools_log_dir = y['source'] + 'aligned/log/'
    filter_bam = y['source'] + 'aligned/' + y['sample'] + '_filtered.bam'
    sorted_bam = y['source'] + 'aligned/' + y['sample'] + '_filtered_sorted.bam'
    clustered_bam = y['source'] + 'aligned/' + y['sample'] + '_clustered.bam'
    clustered_sam = y['source'] + 'aligned/' + y['sample'] + '_clustered.sam'
    stats = y['source'] + 'aligned/' + y['sample'] + '_deduplicated'
    coverage_file = y['source'] + 'aligned/' + y['sample'] + '_coverage.tsv'

    # umi cluster process
    umi_depth_to_cnt = defaultdict(int)
    depth_distribution = [0]

    logger_umitools_process, logger_umitools_errors = store_logs_align(umitools_log_dir)
    umitools_return = umitools(y['samtools_path'], y['umitools_path'], filter_sam, filter_bam, sorted_bam, clustered_bam,
                               clustered_sam, stats, logger_umitools_process, logger_umitools_errors)
    for row in umitools_return.decode('utf-8').split('\n'):
        if not row.startswith('#'):
            if '#umis' in row:
                num_mts = row.split(' ')[-1].strip()
            if 'Input Reads' in row:
                num_input = int(row.split(' ')[-1].strip())
            if 'reads out' in row:
                num_output = int(row.split(' ')[-1].strip())
    num_consilidated = num_input - num_output

    (total_num_bases, mean_bd, min_bd, max_bd, num_50x, num_100x, num_200x, num_500x,
     num_5mean, num_25mean, num_50mean, num_75mean, num_mean) \
        = stats_target_bases(y['samtools_path'], sorted_bam, y['ref_bed'], coverage_file)

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

    print ('Number of total MTs == {0}'.format(num_mts))
    print ('Number of consilidated reads == {0}'.format(num_consilidated))
    print ('Number of total bases on-target == {0}'.format(total_num_bases))
    print ('Mean of base depth == {0}'.format(mean_bd))
    print ('Finish.')

    print('#################################')
    print('7.reformat samfile is starting')
    print('#################################')

    # reformat samfile parameters
    vcready_sam = y['source'] + 'aligned/' + y['sample'] + '_vcready.sam'

    # reformat samfile process
    reformat_sam(clustered_sam, vcready_sam)
    vcready_bam = sam_to_bam(vcready_sam, y['samtools_path'])
    vcready_sorted_bam = sort_index(vcready_bam, y['samtools_path'])
    print ('Finish.')

    print('#################################')
    print('8.variant calling is starting')
    print('#################################')

    # variant calling parameters
    call_outprefix = y['source'] + y['sample'] + '_variant'
    call_log = y['source'] + y['sample'] + '_logfile'

    # variant calling process
    call(y['smcounter_path'], call_outprefix, vcready_sorted_bam, y['normal'], y['ref_bed'], y['ref_fasta'], str(y['rpb']),
         str(y['n_cpu']), str(y['min_bq']), str(y['min_mapq']), str(y['hp_len']), str(y['mismatch']), str(y['mt_drop']),
         str(y['threshold']), str(y['min_frequency']), str(y['min_active_score']), str(y['tlod_threshold']),
         str(y['nlod_threshold']), y['ref_genome'], y['bed_tandem_repeats'], y['bed_repeat_masker_subset'],
         y['bedtools_path'], call_log)
    print ('Finish.')

    print('#################################')
    print('9.annotation is starting')
    print('#################################')

    # annotation parameters
    variant_vcf = y['source'] + y['sample'] + '_variant.smCounter.somatic.cut.vcf'
    annotated_csv = y['source'] + y['sample'] + '_annotated.csv'
    annotation_stats_file = y['source'] + y['sample'] + '_annotate_stats.txt'

    # annotation process
    cosmic_ds, clinvar_ds, g1000_ds = fetch_mysql()
    annotation(cosmic_ds,clinvar_ds,g1000_ds,variant_vcf,annotated_csv,annotation_stats_file)
    fill_table(annotated_csv, y['geneinfo'])
    print ('Finish.')

    print('#################################')
    print('10.statistics is starting')
    print('#################################')

    # statistics parameters
    pipeline_stats_file = y['source'] + y['sample'] + '_statistics_in_report.csv'

    # statistics process
    with open(pipeline_stats_file, 'w') as f:
        f.write('library name,' + y['sample'] + '\n')
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

    print ('Finish.')
    print ('The time of pipeline cost is %s minutes.' % str((time.time() - time_start) / 60))


if __name__ == '__main__':
    main()
