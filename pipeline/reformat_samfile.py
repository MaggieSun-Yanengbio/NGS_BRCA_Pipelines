__author__ = 'Maggie Ruimin Sun'
__version__ = 'v0.2'

import logging
import os
import sys
import time

hg19_chr_length = {
        'chr1': 249250621,
        'chr2': 243199373,
        'chr3': 198022430,
        'chr4': 191154276,
        'chr5': 180915260,
        'chr6': 171115067,
        'chr7': 159138663,
        'chr8': 146364022,
        'chr9': 141213431,
        'chr10': 135534747,
        'chr11': 135006516,
        'chr12': 133851895,
        'chr13': 115169878,
        'chr14': 107349540,
        'chr15': 102531392,
        'chr16': 90354753,
        'chr17': 81195210,
        'chr18': 78077248,
        'chr19': 59128983,
        'chr20': 63025520,
        'chr21': 48129895,
        'chr22': 51304566,
        'chrX': 155270560,
        'chrY': 59373566,
        'chrM': 16571,
    }

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def store_logs_reformat(log_dir):
    formatter_reformat_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_reformat_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_reformat_process = setup_logger('Running Messages of reformating sam file',
                                           log_dir + '/reformat_process.log',
                                           formatter_reformat_process)
    logger_reformat_errors = setup_logger('Errors & Warnings of reformating sam file',
                                          log_dir + '/reformat_errors.log',
                                          formatter_reformat_errors)

    return logger_reformat_process, logger_reformat_errors


def reformat_sam(alignment_sam, output_sam):
    sam = open(alignment_sam)
    sam.readline()
    sam_out = open(output_sam, 'w')
    for chrom in sorted(hg19_chr_length):
        sam_out.write('@SQ\tSN:' + chrom + '\tLN:' + str(hg19_chr_length[chrom]) + '\n')
    for row in sam:
        if row[0:3] == '@SQ':
            continue
        if row[0:3] == '@PG':
            sam_out.write(row)
            continue
        qname, flag, rname, pos, mapq, cigar, rmate, pmate = row.strip().split()[0:8]
        qname, umi = qname.split('_')
        chrom, start = rname.split('_')[0:2]
        pos = str(int(pos) + int(start))
        pmate = str(int(pmate) + int(start))
        if flag == '99' or flag == '163':
            qname += ':' + chrom + '-0-' + pos + '-' + umi + ':' + umi
        elif flag == '83' or flag == '147':
            qname += ':' + chrom + '-1-' + pos + '-' + umi + ':' + umi
        else:

            continue  # error flags
        sam_out.write('\t'.join([qname, flag, chrom, pos, mapq,
                                 cigar, rmate, pmate]) + '\t' + '\t'.join(row.split()[8:]) + '\n')
    sam.close()
    sam_out.close()


def sam_to_bam(samfile,samtools):
    bam = samfile[:-3] + 'bam'
    cmd1 = samtools + ' view -b -S ' + samfile + ' > ' + bam
    os.system(cmd1)

    return bam

def sort_index(bam,samtools):
    sorted = bam[:-4] + '_sorted.bam'
    cmd2 = samtools + ' sort ' + bam + ' > ' + sorted
    os.system(cmd2)
    cmd3 = samtools + ' index ' + sorted
    os.system(cmd3)

    return sorted