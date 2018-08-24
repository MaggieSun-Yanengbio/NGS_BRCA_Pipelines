import os
import sys
from subprocess import Popen, PIPE

def stats_target_bases(samtools_path, filtered_bam, bedfile, cov_file):
    bed = open(bedfile, 'r')
    target_bases_depth = []
    total_num_bases = 0
    num_target_bases = 0
    num_50x = 0
    num_100x = 0
    num_200x = 0
    num_500x = 0
    num_5mean = 0
    num_25mean = 0
    num_50mean = 0
    num_75mean = 0
    num_mean = 0
    bed.readline()
    for row in bed.readlines():
        chr, start, end = row.strip().split('\t')[:3]
        cmd = samtools_path + ' depth ' + filtered_bam + ' -r ' + '_'.join([chr, start, end]) + ':1-' + \
              str(int(end) - int(start)) + ' >> ' + cov_file
        os.system(cmd)
    for row in open(cov_file, 'r'):
        region, pos, cnt = row.strip().split('\t')
        total_num_bases += int(cnt)
        num_target_bases += 1
        target_bases_depth.append(int(cnt))
    mean_bd = round(total_num_bases/num_target_bases, 3)
    min_bd = min(target_bases_depth)
    max_bd = max(target_bases_depth)
    for dp in target_bases_depth:
        if dp >= 50:
            num_50x += 1
            if dp >= 100:
                num_100x += 1
                if dp >= 200:
                    num_200x += 1
                    if dp >= 500:
                        num_500x += 1

        if dp >= mean_bd * 0.05:
            num_5mean += 1
            if dp >= mean_bd * 0.25:
                num_25mean += 1
                if dp >= mean_bd * 0.5:
                    num_50mean += 1
                    if dp >= mean_bd * 0.75:
                        num_75mean += 1
                        if dp >= mean_bd:
                            num_mean += 1

    return (total_num_bases,
            mean_bd,
            min_bd,
            max_bd,
            round(100 * num_50x / num_target_bases, 3),
            round(100 * num_100x / num_target_bases, 3),
            round(100 * num_200x / num_target_bases, 3),
            round(100 * num_500x / num_target_bases, 3),
            round(100 * num_5mean / num_target_bases, 3),
            round(100 * num_25mean / num_target_bases, 3),
            round(100 * num_50mean / num_target_bases, 3),
            round(100 * num_75mean / num_target_bases, 3),
            round(100 * num_mean / num_target_bases, 3))



def umitools(samtools_path, umitools_path, filein, bam, sorted_bam, clustered_bam, clustered_sam, stats):
    cmd1 = samtools_path + ' view -bS ' + filein + ' > ' + bam
    os.system(cmd1)
    cmd2 = samtools_path + ' sort ' + bam + ' > ' + sorted_bam
    os.system(cmd2)
    cmd3 = samtools_path + ' index ' + sorted_bam
    os.system(cmd3)
    process = Popen(['python3', umitools_path, 'dedup', '-I', sorted_bam, '--output-stats=%s'%stats, '-S', clustered_bam,
                     '--edit-distance-threshold', '2', '--paired'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    cmd4 = samtools_path + ' view -h ' + clustered_bam + ' > ' + clustered_sam
    os.system(cmd4)

    return stdout