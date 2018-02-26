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

def store_logs(log_dir):
    formatter_trim_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_trim_errors = logging.Formatter("%(asctime)s;%(levelname)s;                                             %(message)s")
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
			l1 = str(f.readline(),'utf-8')
			if not l1:
				break
			l2 = str(f.readline(),'utf-8')
			l3 = str(f.readline(),'utf-8')
			l4 = str(f.readline(),'utf-8')
			yield [l1, l2, l3, l4]

def trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len, stats_file,
	                logger_trim_process, logger_trim_errors):
	num_total_reads = 0
	num_short_reads = 0
	num_error_reads = 0
	fout1 = open(trimmed1, 'w')
	fout2 = open(trimmed2, 'w')

	for r1, r2 in zip(read_fq(read1,logger_trim_process,logger_trim_errors), 
		read_fq(read2,logger_trim_process,logger_trim_errors)):
		num_total_reads += 1
		if r1[0][0] != '@' or r2[0][0] != '@':
			num_error_reads += 1
			logger_trim_process.info("Error read pair: \n\t"
				+'\t'.join(r1)+'\n\t'+'\t'.join(r2))
		else:
			umi = r2[1][0:12]+';'+r2[3][0:12] # the first 12bp is molecular barcode sequence
			seq2 = r2[1][23:] # the next 11bp is the common sequence
			qual2 = r2[3][23:]
			if len(seq2) < min_read_len:
				num_short_reads += 1
				logger_trim_process.info("Short read pair: \n\t"
				+'\t'.join(r1)+'\n\t'+'\t'.join(r2))
			else:
				h1 = r1[0].strip()+' '+umi+'\n'
				h2 = r2[0].strip()+' '+umi+'\n'
				fout1.write(h1+r1[1]+r1[2]+r1[3])
				fout2.write(h2+seq2+r2[2]+qual2)
	fout1.close()
	fout2.close()
	print('Total number of reads == '+str(num_total_reads))
	print('Number of short reads (<{0}bp) == {1}'.format(min_read_len, num_short_reads))
	print('Number of error reads == '+str(num_error_reads))

	stats_out = open(stats_file, 'w')
	stats_out.write('Total number of reads == '+str(num_total_reads)+'\n')
	stats_out.write('Number of short reads (<{0}bp) == {1}\n'.format(
		min_read_len, num_short_reads))
	stats_out.write('Number of error reads == '+str(num_error_reads)+'\n')

def main():
	if len(sys.argv) < 5:
		sys.exit("Error: Incorrect arguments!")
	source, sample, tailname, out_dir = sys.argv[1:5]
	if len(sys.argv) == 6:
		min_read_len = int(sys.argv[5])
	else:
		min_read_len = 35
	log_dir = out_dir + 'log/'
	if not os.path.exists(log_dir):
		os.makedirs(log_dir)
	tailname = '_'+tailname
	read1 = source + sample + tailname + '_R1.fastq.gz'
	read2 = source + sample + tailname + '_R2.fastq.gz'
	trimmed1 = out_dir + sample + '_R1_undetermined.fastq'
	trimmed2 = out_dir + sample + '_R2_undetermined.fastq'
	stats_file = out_dir + sample + '_basic_stats.txt'
	logger_trim_process, logger_trim_errors = store_logs(log_dir)
	trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len, stats_file,
		logger_trim_process, logger_trim_errors)

if __name__ == '__main__':
	main()