__author__ = 'Maggie Ruimin Sun'
__version__ = 'v0.2'

import os
import sys
import re
import logging
import time
import numpy as np
from difflib import SequenceMatcher

time_start = time.time()
def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def store_logs_filter(log_dir):
    formatter_filter_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_filter_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_filter_process = setup_logger('Post-alignment Filtration Messages',
                                         log_dir + '/filter_process.log',
                                         formatter_filter_process)
    logger_filter_errors = setup_logger('Errors & Warnings of filtration',
                                        log_dir + '/filter_errors.log',
                                        formatter_filter_errors)
    return logger_filter_process, logger_filter_errors


def filter_alignment_samtools(samtools_dir, alignment_sam, min_mapq,
                              max_soft_clip, out_file, stats_file,
                              logger_filter_process, logger_filter_errors):
    stats_file_tmp = stats_file + '.tmp'
    command_count = '{0} view {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, alignment_sam, stats_file_tmp)
    logger_filter_process.info('Samtools counts the total number of read pairs.')
    os.system(command_count)

    # Count supplimentary alignments
    command_count1 = '{0} view -f 2048 {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, alignment_sam, stats_file_tmp)
    logger_filter_process.info('Samtools counts the supplimentary alignments.')
    os.system(command_count1)

    out_file_tmp1 = out_file + '_tmp1.sam'
    command_samtools = '{0} view -ShF 2048 {1} > {2}'.format(
        samtools_dir, alignment_sam, out_file_tmp1)
    os.system(command_samtools)
    # Count unmapped reads
    command_count2 = '{0} view -f 8 {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, out_file_tmp1, stats_file_tmp)
    logger_filter_process.info('Samtools counts the total number of unmapped read pairs.')
    os.system(command_count2)
    # Discard all unmapped read pairs
    out_file_tmp2 = out_file + '_tmp2.sam'
    command_samtools = '{0} view -ShF 8 {1} > {2}'.format(samtools_dir, out_file_tmp1, out_file_tmp2)
    os.system(command_samtools)
    command_samtools = '{0} view -ShF 4 {1} > {2}'.format(samtools_dir, out_file_tmp2, out_file_tmp1)
    os.system(command_samtools)
    # Discard read pairs mapped to different target sequences
    command_samtools1 = "{0} view -Sh {1} \
    | perl -lane 'print if $F[0] =~ /^@/; print if $F[6] =~ /=/;' > {2}".format(
        samtools_dir, out_file_tmp1, out_file_tmp2)
    logger_filter_process.info('Samtools discards the secondary/supplimentary/unmapped alignments\
        and those read pairs mapped to different target regions.')
    os.system(command_samtools1)
    command_count3 = "{0} view {1} | cut -f1 | uniq | wc -l >> {2}".format(
        samtools_dir, out_file_tmp2, stats_file_tmp)
    os.system(command_count3)
    # Discard alignments with MAPQ < min_mapq
    command_count4 = '{0} view -q {1} {2} | cut -f1 | uniq | wc -l >> {3}'.format(
        samtools_dir, min_mapq, out_file_tmp2, stats_file_tmp)
    logger_filter_process.info('Samtools counts the alignments with higher MapQ.')
    os.system(command_count4)
    # Drop the above alignments and those beginning with exact matches > max_soft_clip
    # command_samtools2 = "{0} view -Shq {1} {2}\
    # | perl -lane 'print if $F[0] =~ /^@/; next unless $F[5] =~ /^(\d+)M/;print if $1 >= {3}'\
    # >{4}".format(samtools_dir, min_mapq, out_file_tmp2, max_soft_clip, out_file)
    command_samtools2 = "{0} view -Shq {1} {2} > {3}".format(samtools_dir, min_mapq,
                                                             out_file_tmp2, out_file)
    os.system(command_samtools2)
    logger_filter_process.info('Samtools discards the alignments with lower MapQ and mismatched beginning.')
    command_count5 = '{0} view {1} | cut -f1 | uniq | wc -l >> {2}'.format(samtools_dir, out_file, stats_file_tmp)
    os.system(command_count5)
    print('Compeleted filtering alignments with samtools.')

    # Output the alignment statistics.
    stats = open(stats_file_tmp)
    (total_count, sec_count, unmap_count, same_chr_count,
     hmapq_count, final_count) = [int(x.strip()) for x in stats.readlines()[0:6]]
    stats.close()
    stats_out = open(stats_file, 'w')
    stats_out.write('Total number of read pairs == {0}\n'.format(total_count))
    stats_out.write('Number of secondary/supplimentary alignments == {0}\n'.format(sec_count))
    stats_out.write('Number of unmapped read pairs == {0}\n'.format(unmap_count))
    remained = total_count - sec_count - unmap_count
    stats_out.write('Number of read pairs mapped to different target regions == {0}\n'.format(
        remained - same_chr_count))
    stats_out.write('Number of low MAPQ alignments (< {0}) == {1}\n'.format(min_mapq, same_chr_count - hmapq_count))
    # stats_out.write('Number of alignments beginning with less than {0} exact matches == {1}\n'.format(
    #    max_soft_clip, hmapq_count-final_count))
    stats_out.write('Number of alignments passed SAMTools filtration == {0}\n'.format(final_count))
    stats_out.write(
        'Percent of alignments passed SAMTools filtration == {0}%\n'.format(100 * final_count / total_count))
    stats_out.write('Time cost at samtools filtration == {0} min\n'.format((time.time() - time_start) / 60))
    stats_out.close()
    os.system('rm ' + stats_file_tmp)
    os.system('rm ' + out_file_tmp1 + ' ' + out_file_tmp2)
    return (sec_count, unmap_count, min_mapq, same_chr_count - hmapq_count, final_count, round(100 * final_count / total_count, 3))


def complement_base(base):
    if base == 'A':
        return 'T'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'T':
        return 'A'
    else:
        return base


def reverse_complement(seq):
    seq_rc = ''
    ls = len(seq)
    for i in range(ls - 1, -1, -1):
        seq_rc += complement_base(seq[i])
    return seq_rc


def compare_seqs(primer, seq, strand):
    matcher = SequenceMatcher(None, primer, seq)
    dist = len(primer)
    i = 0
    j = 0
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'equal':
            if i == 0:
                dist -= i1
                if strand == 0:
                    match_start = j1
            dist -= (i2-i1)
            j = j2
            i += 1
    if strand == 1:
        match_start = len(seq) - j
    return dist, match_start


def identify_gs_primers(samtools_dir, alignment_sam, primers_file, max_dist,
                        out_file, stats_file, primer_stats_file,
                        logger_filter_process, logger_filter_errors):
    primers = {}
    primer_pos = {}
    csv = open(primers_file)
    # csv.readline()
    for row in csv:
        chrom, strand, seq, start, stop, gene = row.strip().split(',')[0:6]
        if chrom == 'chr':
            continue
        pos_start = int(start)
        pos_stop = int(stop)
        if chrom not in primer_pos:
            primer_pos[chrom] = {}

        if strand == '0':
            for pos in range(pos_start - 100, pos_stop + 101):
                if pos not in primer_pos[chrom]:
                    primer_pos[chrom][pos] = []
                primer_pos[chrom][pos].append((seq, chrom + '_' + start + '_' + stop))
        else:
            for pos in range(pos_stop - 100, pos_start + 101):
                if pos not in primer_pos[chrom]:
                    primer_pos[chrom][pos] = []
                primer_pos[chrom][pos].append((reverse_complement(seq), chrom + '_' + stop + '_' + start))

        if strand == '0':
            key_primer = chrom + '_' + start + '_' + stop
            primers[key_primer] = {}
            primers[key_primer]['strand'] = '0'
            primers[key_primer]['seq'] = seq
        else:
            key_primer = chrom + '_' + stop + '_' + start
            primers[key_primer] = {}
            primers[key_primer]['strand'] = '1'
            primers[key_primer]['seq'] = reverse_complement(seq)
        primers[key_primer]['gene'] = gene
        primers[key_primer]['chromosome'] = chrom
        primers[key_primer]['start'] = start
        primers[key_primer]['stop'] = stop
        primers[key_primer]['length'] = len(seq)
        primers[key_primer]['on-target-reads'] = 0
        primers[key_primer]['on-target-mts'] = set()
        # primers[key_primer]['mis-target'] = []

    csv.close()
    num_primers = len(primers)
    print('Number of primers of interest == ' + str(num_primers))

    num_total_alignments = 0
    num_short_fragment = 0
    num_on_target = 0
    num_off_target = 0
    num_off_target_wrong_strand = 0
    num_unproper_pairs = 0
    num_target_primers = 0
    num_primers_5mean_rdp = 0
    num_primers_25mean_rdp = 0
    num_primers_50mean_rdp = 0
    num_primers_75mean_rdp = 0
    num_primers_mean_rdp = 0
    num_primers_5mean_mdp = 0
    num_primers_25mean_mdp = 0
    num_primers_50mean_mdp = 0
    num_primers_75mean_mdp = 0
    num_primers_mean_mdp = 0
    primer_to_mts = {}
    off_target_reads = {}

    soft_clips_start = re.compile('^(\d+)S')
    soft_clips_end = re.compile('(\d+)S$')
    flag_eligible = ['83', '99']

    sam = open(alignment_sam)
    ready = open(out_file, 'w')
    off_sam = open(alignment_sam + '.off', 'w')
    proper_pair = True
    sam_rows = sam.readlines()
    while True:
        if sam_rows[0][0] == '@':
            ready.write(sam_rows[0])
            sam_rows.remove(sam_rows[0])
        else:
            break

    sam_rows = sorted(sam_rows)
    i = 0
    while i < len(sam_rows) - 1:
        row = sam_rows[i]
        mate = sam_rows[i + 1]
        qname, flag, rname, pos, mapq, cigar, rmn, pmn, insl, seq = row.strip().split()[0:10]
        umi = qname[-12:]
        if qname in mate:
            proper_pair = True
            i += 2
        else:
            proper_pair = False
            i += 1

        if flag not in flag_eligible:
            if int(flag) > 128:
                flag_mate = mate.split()[1]
                if flag_mate in flag_eligible:
                    qname, flag, rname, pos, mapq, cigar, rmn, pmn, insl, seq = mate.strip().split()[0:10]
                    row, mate = mate, row
            else:
                proper_pair = False

        if not proper_pair:
            num_unproper_pairs += 1
            continue

        num_total_alignments += 1
        chrom, start, stop = rname.split('_')
        if chrom not in primer_pos:
            num_off_target += 1
            continue
        pos = int(pos) + int(start) - 1
        pos_rv = pos + len(seq) - 1

        if soft_clips_start.search(cigar):
            pos_rv = pos_rv - int(soft_clips_start.search(cigar).groups()[0])
        if soft_clips_end.search(cigar):
            pos_rv = pos_rv - int(soft_clips_end.search(cigar).groups()[0])

        off_target = 0
        off = False
        min_dist = max_dist
        match_start0 = len(seq)
        if flag == '83':  # read 1 is mapped to the reversed strand
            if pos_rv not in primer_pos[chrom]:
                off = True
                num_off_target += 1
                num_off_target_wrong_strand += 1
                continue
            for primer_info in primer_pos[chrom][pos_rv]:
                primer_name = primer_info[1]
                primer_seq = primer_info[0]
                lpr = primers[primer_name]['length']
                if primers[primer_name]['strand'] == '0':
                    off_target += 1
                    continue
                dist, match_start = compare_seqs(primer_seq, seq, 1)
                if dist > max_dist:
                    off_target += 1
                elif dist <= min_dist and match_start < match_start0:
                    min_dist = dist
                    match_start0 = match_start
                    primer_on_target = primer_name

            if off_target == len(primer_pos[chrom][pos_rv]):
                off = True
        elif flag == '99':
            if pos in primer_pos[chrom]:
                for primer_info in primer_pos[chrom][pos]:
                    primer_name = primer_info[1]
                    primer_seq = primer_info[0]
                    lpr = primers[primer_name]['length']
                    if primers[primer_name]['strand'] == '1':
                        off_target += 1
                        continue
                    dist, match_start = compare_seqs(primer_seq, seq, 0)
                    if dist > max_dist:
                        off_target += 1
                    elif dist <= min_dist and match_start < match_start0:
                        min_dist = dist
                        match_start0 = match_start
                        primer_on_target = primer_name

                if off_target == len(primer_pos[chrom][pos]):
                    off = True
            else:
                num_off_target += 1
                num_off_target_wrong_strand += 1
                continue
        if off:
            num_off_target += 1
            off_target_reads[rname + '_' + str(pos)] = seq
            off_sam.write(row)
            off_sam.write(mate)
        else:
            num_on_target += 1
            primers[primer_on_target]['on-target-reads'] += 1
            primers[primer_on_target]['on-target-mts'].add(umi)
            ready.write(row)
            ready.write(mate)
    sam.close()
    ready.close()

    ratio_off = round(100 * num_off_target / num_total_alignments, 3)
    ratio_on = round(100 * num_on_target / num_total_alignments, 3)
    print('Total number of alignments (PE) == ' + str(num_total_alignments))
    print('Number of reads mapped in unproper pairs == ' + str(num_unproper_pairs))
    print('Number of off-target alignments == {0} ({1}%)({2})'.format(num_off_target, ratio_off, num_off_target_wrong_strand))
    print('Number of on-target alignments == {0} ({1}%)'.format(num_on_target, ratio_on))

    stats_out = open(stats_file, 'a')
    stats_out.write('Number of unproperly-paired alignments == ' + str(num_unproper_pairs) + '\n')
    stats_out.write('Number of panel primers == ' + str(num_primers) + '\n')
    stats_out.write('Number of off-target alignments == {0} ({1}%)\n'.format(num_off_target, ratio_off))
    stats_out.write('Number of on-target alignments == {0} ({1}%)\n'.format(num_on_target, ratio_on))
    stats_out.close()

    st_out = open(primer_stats_file, 'w')
    st_out.write('chrm,F/R,start,stop,sequence,#on-target-reads,%on-target-reads\n')
    mean_rdp = round(num_on_target / num_primers, 3)
    for pr in sorted(primers):
        on_target_reads = primers[pr]['on-target-reads']
        if on_target_reads > 0:
            num_target_primers += 1
            if on_target_reads > mean_rdp * 0.05:
                num_primers_5mean_rdp += 1
                if on_target_reads >= mean_rdp * 0.25:
                    num_primers_25mean_rdp += 1
                    if on_target_reads >= mean_rdp * 0.5:
                        num_primers_50mean_rdp += 1
                        if on_target_reads >= mean_rdp * 0.75:
                            num_primers_75mean_rdp += 1
                            if on_target_reads >= mean_rdp:
                                num_primers_mean_rdp += 1
        primer_to_mts[pr] = len(primers[pr]['on-target-mts'])
        ratio = 100 * primers[pr]['on-target-reads'] / float(num_on_target)
        st_out.write(','.join([primers[pr]['chromosome'], primers[pr]['strand'],
                               primers[pr]['start'], primers[pr]['stop'],
                               primers[pr]['seq'], str(primers[pr]['on-target-reads']),
                               str(ratio)]) + '\n')
    mean_mdp = round(np.mean(list(primer_to_mts.values())), 3)
    for mts in primer_to_mts.values():
        if mts > mean_mdp * 0.05:
            num_primers_5mean_mdp += 1
            if mts > mean_mdp * 0.25:
                num_primers_25mean_mdp += 1
                if mts > mean_mdp * 0.5:
                    num_primers_50mean_mdp += 1
                    if mts > mean_mdp * 0.75:
                        num_primers_75mean_mdp += 1
                        if mts > mean_mdp:
                            num_primers_mean_mdp += 1

    st_out.close()

    return (num_unproper_pairs, num_primers, num_off_target, ratio_off, num_on_target, ratio_on, num_target_primers,
            mean_rdp, num_primers_5mean_rdp, num_primers_25mean_rdp, num_primers_50mean_rdp, num_primers_75mean_rdp,
            num_primers_mean_rdp, mean_mdp, num_primers_5mean_mdp, num_primers_25mean_mdp, num_primers_50mean_mdp,
            num_primers_75mean_mdp, num_primers_mean_mdp)
