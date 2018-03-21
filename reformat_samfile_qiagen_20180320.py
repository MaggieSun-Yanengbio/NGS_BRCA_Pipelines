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
    'chrY': 59034049,
    'chrM': 16569,
}

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    
    return logger

def store_logs(log_dir):
    formatter_reformat_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_reformat_errors = logging.Formatter("%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_reformat_process = setup_logger('Running Messages of reformating sam file', 
                                      log_dir + '/reformat_process.log', 
                                      formatter_reformat_process)
    logger_reformat_errors = setup_logger('Errors & Warnings of reformating sam file', 
                                    log_dir + '/reformat_errors.log',
                                    formatter_reformat_errors)
    return logger_reformat_process, logger_reformat_errors

"""
The numbers of consolidated sequences in R1 and R2 are different, because those
consolidated sequences with too many N's are removed. 
"""
def filter_consolidation(consolidated_read1, consolidated_read2, sample_id_list):
    consolidate_dict = {}
    fq1 = open(consolidated_read1)
    reads1 = fq1.readlines()
    fq1.close()

    fq2 = open(consolidated_read2)
    reads2 = fq2.readlines()
    fq2.close()
    num_reads1 = len(reads1)
    num_reads2 = len(reads2)
    print([num_reads1/4, num_reads2/4])

    temp_dict = {}
    i = 0
    while i < num_reads1:
        header = reads1[i]
        seq = reads1[i+1]
        qual = reads1[i+3]
        i += 4
        sample_id, umi = header.split(';')[0:2]
        temp_dict[umi] = {}
        temp_dict[umi]['cons_sample_id'] = sample_id
        temp_dict[umi]['R1'] = [header, seq, '+\n', qual]

    j = 0
    while j < num_reads2:
        header = reads2[j]
        sample_id, umi = header.split(';')[0:2]
        if umi in temp_dict:
            temp_dict[umi]['R2'] = [header, reads2[j+1], reads2[j+2], reads2[j+3]]
        j+= 4

    for key_id in list(temp_dict):
        if len(temp_dict[key_id]) < 3:
            rm_key = temp_dict.pop(key_id, None)

    print('Number of overlapping consolidated read pairs: '+str(len(temp_dict)))
    
    fout1 = open(consolidated_read1+'_filtered.fastq', 'w')
    fout2 = open(consolidated_read2+'_filtered.fastq', 'w')
    for key_id in temp_dict:
        fout1.write(''.join(temp_dict[key_id]['R1']))
        fout2.write(''.join(temp_dict[key_id]['R2']))
    fout1.close()
    fout2.close()

    fs = open(sample_id_list)
    for line in fs:
        sample_id, umi = line.strip('\t').split()
        sample_id = sample_id[1:]
        if umi in temp_dict:
            consolidate_dict[sample_id] = {}
            consolidate_dict[sample_id]['umi'] = umi
            consolidate_dict[sample_id]['cons_sample_id'] = temp_dict[umi]['cons_sample_id']
    fs.close()
    return consolidate_dict

def reformat_sam(consolidated_read1, consolidated_read2, sample_id_list, 
    alignment_sam, output_sam):
    consolidate_dict = filter_consolidation(consolidated_read1, consolidated_read2, sample_id_list)
    sam = open(alignment_sam)
    sam_out = open(output_sam, 'w')
    for chrom in hg19_chr_length:
        sam_out.write('@SQ\tSN:'+chrom+'\tLN:'+str(hg19_chr_length[chrom])+'\n')
    for row in sam:
        if row[0:3] == '@SQ':
            continue
        if row[0:3] == '@PG':
            sam_out.write(row)
            continue
        qname, flag, rname, pos, mapq, cigar, rmate, pmate = row.strip().split()[0:8]
        if qname not in consolidate_dict:
            continue
        umi = consolidate_dict[qname]['umi']
        chrom, start = rname.split('_')[0:2]
        pos = str(int(pos) + int(start) - 1)
        pmate = str(int(pmate) + int(start) - 1)
        if flag == '99' or flag == '163':
            qname += ':'+chrom+'-0-'+pos+'-'+umi+':'+umi
        elif flag == '83' or flag == '147':
            qname += ':'+chrom+'-1-'+pos+'-'+umi+':'+umi
        else:

            continue #error flags
        sam_out.write('\t'.join([qname, flag, chrom, pos, mapq, 
            cigar, rmate, pmate])+'\t'+'\t'.join(row.split()[8:])+'\n')
    sam.close()
    sam_out.close()

def main():
    if len(sys.argv) < 3:
        print('python3 reformat_samfile_qiagen_20180320.py \
            /path/to/input/ sample_name')
        sys.exit("Error: Incorrect arguments!")
    source, sample_name = sys.argv[1:3]
    consolidated_read1 = source + 'consolidated/'+sample_name+'_R1_consolidated.fastq'
    consolidated_read2 = source + 'consolidated/'+sample_name+'_R2_consolidated.fastq'
    sample_id_list = source + 'undetermined/'+sample_name+'_R1_undetermined.fastq_umi.fastq_ids.txt'
    alignment_sam = source + 'aligned/'+sample_name+'_filtered.sam'
    output_sam = source + 'aligned/'+sample_name+'_vcready.sam'
    
    time_start = time.time()
    reformat_sam(consolidated_read1, consolidated_read2, sample_id_list, alignment_sam, output_sam)
    print('Finish reformating alignment SAM file after {0} min.'.format((time.time()-time_start)/60)) 
if __name__ == '__main__':
    main()
