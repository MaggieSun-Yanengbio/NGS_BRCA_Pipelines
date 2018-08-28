__author__ = 'Maggie Ruimin Sun'
__version__ = 'v0.2'

import logging
import os
import re
import sys
import time
import gzip
import itertools

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def store_logs_trim(log_dir):
    formatter_trim_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_trim_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_trim_process = setup_logger('Running Messages of read trimming',
                                       log_dir + '/trim_process.log',
                                       formatter_trim_process)
    logger_trim_errors = setup_logger('Errors & Warnings of read trimming',
                                      log_dir + '/trim_errors.log',
                                      formatter_trim_errors)
    return logger_trim_process, logger_trim_errors


def read_fq(file_name, logger_trim_process, logger_trim_errors):
    if not os.path.isfile(file_name):
        logger_trim_errors.error("%s does not exist!\n", file_name)
        print(file_name + ' does not exist!')
    if re.search('.gz$', file_name):
        fastq = gzip.open(file_name, 'r')
    else:
        fastq = open(file_name)

    with fastq as f:
        while True:
            l1 = str(f.readline(), 'utf-8')
            if not l1:
                break
            l2 = str(f.readline(), 'utf-8')
            l3 = str(f.readline(), 'utf-8')
            l4 = str(f.readline(), 'utf-8')
            yield [l1, l2, l3, l4]


def trim_read1(r1, common_seq2, mt_barcode):
    l1 = len(r1[1].strip())
    r1_end3 = r1[1].strip()[(l1 - 23):]
    trim_seq = common_seq2 + mt_barcode
    for i in range(23):
        if r1_end3[i:] == trim_seq[0:(23 - i)]:
            break
    pos_trim_r1 = l1 - 23 + i
    return pos_trim_r1, [r1[0], r1[1].strip()[0:pos_trim_r1] + '\n', r1[2], r1[3].strip()[0:pos_trim_r1] + '\n']


def trim_read2(r2, common_seq1):
    l2 = len(r2[1].strip())
    r2_end3 = r2[1].strip()[(l2 - 21):]
    for i in range(21):
        if r2_end3[i:] == common_seq1[0:(21 - i)]:
            break
    pos_trim_r2 = l2 - 21 + i
    return pos_trim_r2, [r2[0], r2[1].strip()[0:pos_trim_r2] + '\n', r2[2], r2[3].strip()[0:pos_trim_r2] + '\n']


def trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len, common_seq1,
                    common_seq2, stats_file, logger_trim_process, logger_trim_errors):
    time_start = time.time()
    num_total_reads = 0
    num_short_reads = 0
    num_error_reads = 0
    num_unproper_reads = 0
    num_without_comseq_reads = 0
    num_clean_reads = 0
    fout1 = open(trimmed1, 'w')
    fout2 = open(trimmed2, 'w')

    for r1, r2 in zip(read_fq(read1, logger_trim_process, logger_trim_errors),
                      read_fq(read2, logger_trim_process, logger_trim_errors)):
        num_total_reads += 1
        if r1[0][0] != '@' or r2[0][0] != '@':
            num_error_reads += 1
            num_unproper_reads += 1
            logger_trim_process.info("Error read pair: \n\t"
                                     + '\t'.join(r1) + '\n\t' + '\t'.join(r2) + '\n')
        else:
            start_common = r2[1].find(common_seq2)
            if start_common < 12:
                num_error_reads += 1
                num_without_comseq_reads += 1
                logger_trim_process.info("Error barcode/common seqs:" + str(start_common) + "\n\t"
                                         + '\t'.join(r1) + '\n\t' + '\t'.join(r2) + '\n')
            else:
                umi = r2[1][(start_common - 12):start_common]
                qua = r2[3][(start_common - 12):start_common]
                r2[1] = r2[1][(start_common + 11):]
                r2[3] = r2[3][(start_common + 11):]
                pos_trim_r1, r1 = trim_read1(r1, common_seq2, umi)
                pos_trim_r2, r2 = trim_read2(r2, common_seq1)
                if pos_trim_r1 < min_read_len or pos_trim_r2 < min_read_len:
                    num_short_reads += 1
                    logger_trim_process.info("Short read pair: \n\t"
                                             + '\t'.join(r1) + '\n\t' + '\t'.join(r2) + '\n')
                else:
                    num_clean_reads += 1
                    h1 = r1[0].split(' ')[0] + '_' + umi + ' ' + r1[0].split(' ')[1]
                    h2 = r2[0].split(' ')[0] + '_' + umi + ' ' + r2[0].split(' ')[1]
                    fout1.write(h1 + r1[1] + r1[2] + r1[3])
                    fout2.write(h2 + r2[1] + r2[2] + r2[3])
    fout1.close()
    fout2.close()
    print('Total number of reads == ' + str(num_total_reads))
    print('Number of short reads (<{0}bp) == {1}'.format(min_read_len, num_short_reads))
    print('Number of error reads == ' + str(num_error_reads))

    stats_out = open(stats_file, 'w')
    stats_out.write('Total number of reads == ' + str(num_total_reads) + '\n')
    stats_out.write('Number of short reads (<{0}bp) == {1}\n'.format(
        min_read_len, num_short_reads))
    stats_out.write('Number of error reads == ' + str(num_error_reads) + '\n')
    stats_out.write('The time of trimming is %s minutes.' % str((time.time() - time_start) / 60))

    return (num_total_reads, num_short_reads, num_unproper_reads, num_without_comseq_reads,
            num_clean_reads, round(100*num_clean_reads/num_total_reads, 3))

