import os
import re
from itertools import groupby
import shlex
import subprocess
import logging

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def store_logs_qc(log_dir):
    formatter_qc_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_qc_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_qc_process = setup_logger('Running Messages of static',
                                       log_dir + '/static_process.log',
                                       formatter_qc_process)
    logger_qc_errors = setup_logger('Errors & Warnings of pipeline',
                                      log_dir + '/static_errors.log',
                                      formatter_qc_errors)
    return logger_qc_process, logger_qc_errors


def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr

def split_n_bases(nbases):
    if nbases == 'NULL':
        return ["NULL" ,"NULL"]
    else:
        bases = nbases.split(",")
        bases_posi =[]
        for posi in range(0,len(bases)):
            posi = bases[posi]
            if "-" in str(posi):
                for i in range(int(posi.split("-")[0]), int(posi.split("-")[1])+1):
                    bases_posi.append(i)
            else:
                bases_posi.append(int(posi))
        l = bases_posi
        a = list()
        b = []
        result = []
        function = lambda x: x[1] - x[0]
        for k, g in groupby(enumerate(l), function):
            g = list(g)
            b.append(k)
            a.append(g)
        for c in range(len(b)):
            if [v for i, v in a[c]][0] != [v for i, v in a[c]][-1]:
                result.append("%d-%d" % ([v for i, v in a[c]][0], [v for i, v in a[c]][-1]))
            if [v for i, v in a[c]][0] == [v for i, v in a[c]][-1]:
                result.append([v for i, v in a[c]][0])
        if "-" not in str(result[len(result)-1]):
            result_last = result[len(result)-1]
        else:
            result_last = result[len(result)-1].split('-')[1]
    return ["[" + str(result[0]) + "]" , str(result_last)]



def getinfo(fastqc, logger_statistics_process, logger_statistics_errors):
    qc_read = fastqc.split(".fastq")[0] + '_fastqc.zip'
    if not os.path.isfile( fastqc.split(".fastq")[0] + '_fastqc'):
        command1 = 'unzip' + ' -o ' + qc_read  + ' -d ' + os.path.dirname(qc_read)
        stdout, stderr = stdout_err(command1)
        logger_statistics_process.info(stdout)
        logger_statistics_errors.info(stderr)
    # print('{0} has been completed.'.format(qc_read))
    qc_data = qc_read.rstrip(".zip") + '/' + 'fastqc_data.txt'
    qcdata = open(qc_data,"r")
    modules = 0
    # 定义正则表达式
    value = re.compile(r'\d+')
    perbasesequencequalit_posi = []
    perbasesequencequalit = []
    persequencequalityscores = []
    persequencequalityreads = []
    per_basen1 = []
    perbase_ncontent = []
    f = qcdata.readlines()
    for lines in range(0 ,len(f)):
        line = f[lines]
        line = line.strip()
        if '>>END_MODULE' in line:
            modules = modules + 1
        if 'Total Sequences' in line:
            raw_reads = re.findall('\d+',line)
        if 'Sequence length' in line:
            #mode = re.compile(r'\d+-\d+')
            #seq_length = mode.findall(line)
            seq_length = line.split('\t')[1]
            if '-' in seq_length:
                seq_length_min,seq_length_max = seq_length.split('-')
            else:
                seq_length_min = str(seq_length)
                seq_length_max = str(seq_length)
            #print(seq_length)
        if '%GC' in line:
            gc = re.findall('\d+',line)
        if modules == 1 and value.match(line[0]):
            perbasesequencequalit_posi.append(line.split('\t')[0])
            perbasesequencequalit.append(float(line.split('\t')[1]))
        if modules == 3 and value.match(line[0]):
            persequencequalityscores.append(line.split('\t')[0])
            persequencequalityreads.append(float(line.split('\t')[1]))
        if modules == 6 and value.match(line[0]):
            per_basen1.append(line.split('\t')[0])
            perbase_ncontent.append(line.split('\t')[1])
    #--Per base sequence quality
    perbasesequencequalit_mean = sum(perbasesequencequalit)/len(perbasesequencequalit)
    lowqualit_bases= []
    for i in range(0,len(perbasesequencequalit)):
        if float(perbasesequencequalit[i]) < perbasesequencequalit_mean-2:
            lowqualit_bases.append(perbasesequencequalit_posi[i])
    if len(lowqualit_bases) == 0:
        lowqualit_bases.append('NULL')
    #--per sequence quality
    persequencequalityreads1 = []
    for i in range(0,len(persequencequalityreads)):
        persequencequalityreads1.append(sum(persequencequalityreads[i:])/int(raw_reads[0]))
    for i in range(0,len(persequencequalityscores)):
        if persequencequalityscores[i] == '20':
            q20 = persequencequalityreads1[i]
        if persequencequalityscores[i] == '30':
            q30 = persequencequalityreads1[i]
    if '20' not in persequencequalityscores:
        q20 = 1.0000
    #--N content
    nbases = []
    for i in range(0,len(perbase_ncontent)):
        if float(perbase_ncontent[i]) > 0:
            nbases.append(per_basen1[i])
    if len(nbases) == 0:
        nbases.append('NULL')
    #---
    #print(raw_reads)
    #print(seq_length)
    #print(perbasesequencequalit_mean)
    #print([raw_reads[0], seq_length[0], GC[0], str(round(Perbasesequencequalit_mean,3)), ','.join(lowqualitBases), str('%.3f%%' % (Q20 * 100)), str('%.3f%%' % (Q30 * 100)), ','.join(nbases)]
    return [raw_reads[0], seq_length_min, seq_length_max, gc[0], str(round(perbasesequencequalit_mean,3)),
            split_n_bases(','.join(lowqualit_bases))[0], str('%.3f%%' % (q20 * 100)), str('%.3f%%' % (q30 * 100)),
            split_n_bases(','.join(nbases))[0], split_n_bases(','.join(nbases))[1]]


def qc_raw_reads(fastQC_dir, out_dir, sample, read1, read2, logger_statistics_process, logger_statistics_errors):
    qc_read = out_dir + '/' + os.path.basename(read1).split(".fastq")[0] + '_fastqc.zip'

    command1 = '{0} {1} {2} -o {3}'.format(fastQC_dir, read1, read2, out_dir)
    stdout, stderr = stdout_err(command1)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)

    qc_statistics = out_dir + '/' + sample + '_qc.statistics.txt'
    qc_result1 = getinfo(out_dir + '/' + os.path.basename(read1), logger_statistics_process, logger_statistics_errors)
    qc_result2 = getinfo(out_dir + '/' + os.path.basename(read2), logger_statistics_process, logger_statistics_errors)
    min_length1, max_length1, gc1 = qc_result1[1:4]
    q20_r1, q30_r1 = qc_result1[6:8]
    min_length2, max_length2, gc2 = qc_result2[1:4]
    q20_r2, q30_r2 = qc_result2[6:8]
    fout = open(qc_statistics, 'w')
    fout.write('\t'.join(['SampleID', 'Sequence direction', 'raw reads', 'min length', 'max length', 'GC content',
                          'mean of Per base qualit', 'low qualit Bases position', 'Q20', 'Q30',
                          'N_bases position']) + '\n')
    fout.write('\t'.join([sample, read1.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result1[0:9])]) + '\n')
    fout.write('\t'.join([sample, read2.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result2[0:9])]) + '\n')
    fout.close()

    return (min_length1, max_length1, gc1, q20_r1, q30_r1), (min_length2, max_length2, gc2, q20_r2, q30_r2)
