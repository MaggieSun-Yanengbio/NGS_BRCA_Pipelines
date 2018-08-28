__author__ = 'Maggie Ruimin Sun'
__version__ = 'v0.2'

import os
import logging
import time
import sys
from subprocess import Popen, PIPE


def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def store_logs_align(log_dir):
    formatter_bwa_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_bwa_errors = logging.Formatter("%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_bwa_process = setup_logger('BWA Running Messages', 
                                      log_dir + '/bwa_process.log', 
                                      formatter_bwa_process)
    logger_bwa_errors = setup_logger('Errors & Warnings of BWA', 
                                    log_dir + '/bwa_errors.log',
                                    formatter_bwa_errors)
    return logger_bwa_process, logger_bwa_errors

def align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, read1, read2, 
	                out_file, num_threads, logger_bwa_process, 
	                logger_bwa_errors): 
    if not os.path.isfile(read1):
    	logger_bwa_errors.error('%s does not exist!', read1)
    	print("Error: cannot find NGS read file!")

    index_file_extensions = ['.pac', '.amb','.ann','.bwt', '.sa']
    genome_indexed = True
    for extension in index_file_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        idx_process = Popen([bwa_dir, 'index', '-p', ref_index_name, ref_fa_file], stdout=PIPE, stderr=PIPE)
        idx_stdout, idx_stderr = idx_process.communicate()
        idx_process.wait()
        logger_bwa_process.info(idx_stdout)
        logger_bwa_errors.info(idx_stderr)
        print('BWA genome index files have been built.')
    else:
        print('BWA genome index files exist.')

    bwa_align_command = '{0} mem -t {1} {2} {3} {4} > {5}'.format(
        bwa_dir, num_threads, ref_index_name, read1, read2, out_file)
    align_process = Popen(bwa_align_command, shell=True, stdout=PIPE, stderr=PIPE)
    align_process.wait()
    align_stdout, align_stderr = align_process.communicate()
    if align_stdout is not b'':
        logger_bwa_process.info(align_stdout)
    if align_stderr is not b'':
        logger_bwa_errors.info(align_stderr)
    print('BWA alignment has been completed.')
    logger_bwa_process.info('BWA alignment has been completed.')
