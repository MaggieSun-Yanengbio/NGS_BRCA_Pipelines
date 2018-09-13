# smCounter: barcode-aware variant caller
# Chang Xu. 23May2016; online version of 10APR2017
import os
import datetime
import subprocess
import math
import operator
import argparse
import random
import multiprocessing
import traceback
from collections import defaultdict
from collections import Counter

# 3rd party modules
import pysam
import numpy as np
import scipy.stats

# global contants (note that multiprocessing imports this .py file, so do not do much work outside functions)
pcr_error = 1e-6
pcr_no_error = 1.0 - 3e-5
atgc = ('A', 'T', 'G', 'C')


# -------------------------------------------------------------------------------------
# function to calculate posterior probability for each barcode.
# -------------------------------------------------------------------------------------
def cal_prob(oneBC, ref, alt, mutation_freq):
    p0_right_bc = 1.0
    p1_right_bc = 1.0
    for readid, info in oneBC.items():
        base, p_error, pair_order = info[0]
        if base == ref:
            p0_right = 1-p_error
            p1_right = (1-mutation_freq)*(1-p_error) + mutation_freq*p_error/3
        elif base == alt:
            p0_right = p_error/3
            p1_right = mutation_freq*(1-p_error) + (1-mutation_freq)*p_error/3
        else:
            p0_right = p_error/3
            p1_right = p_error/3
        p0_right_bc *= p0_right
        p1_right_bc *= p1_right

    return (p0_right_bc, p1_right_bc, base)


# -------------------------------------------------------------------------------------
# classification somatic and germline variants.
# -------------------------------------------------------------------------------------
def variants_classification(chrom, pos, ref, alt, normal_bamfile, min_bq, refgenome):
    input = pysam.AlignmentFile(normal_bamfile, 'rb')
    refseq = pysam.FastaFile(refgenome)
    divide = 1.0
    alt_cnt = 0
    depth_cnt = 0
    allele_cnt = defaultdict(int)
    is_exist_in_normal = False
    is_triple_alleles = False
    is_zero_coverage_in_normal = False

    if alt.startswith('DEL'):
        ref = ref[0]

    for reads in input.pileup(region = chrom + ':' + pos + ':' + pos, truncate=True, stepper='nofilter'):
        for read in reads.pileups:
            #get base and base quality
            if read.indel > 0:
                site = read.alignment.query_sequence[read.query_position]
                inserted = read.alignment.query_sequence[(read.query_position + 1): (read.query_position + 1 + read.indel)]
                base = 'INS|' + site + '|' + site + inserted
                bq = read.alignment.query_qualities[read.query_position]
            elif read.indel < 0:
                site = read.alignment.query_sequence[read.query_position]
                deleted = refseq.fetch(reference=chrom, start=int(pos), end=int(pos) + abs(read.indel))
                deleted = deleted.upper()
                base = 'DEL|' + site + deleted + '|' + site
                bq = read.alignment.query_qualities[read.query_position]
            else:
                if read.is_del:
                    base = 'DEL'
                    bq = min_bq
                else:
                    base = read.alignment.query_sequence[read.query_position]
                    bq = read.alignment.query_qualities[read.query_position]
            allele_cnt[base] += 1
            p_error = 10**(-float(bq)/10)

            #calculation L(M0)/L(M(m,k))
            if base == ref:
                p0 = 1-p_error
                p_germline = 0.5*(1-p_error) + 0.5*p_error/3
            elif base == alt:
                p0 = p_error / 3
                p_germline = 0.5 * (1 - p_error) + 0.5 * p_error / 3
                alt_cnt += 1
            else:
                p0 = p_error/3
                p_germline = p_error/3
            divide *= p0/p_germline
            depth_cnt += 1
    try:
        nlod = math.log10(divide)
    except ValueError:
        nlod = 0

    if depth_cnt == 0:
        is_zero_coverage_in_normal = True
    else:
        if alt_cnt >= 3 or float(alt_cnt)/depth_cnt >= 0.03:
            is_exist_in_normal = True

        if len(allele_cnt) != 1:
            most_common_allele, second_common_allele = sorted(allele_cnt.items(), key=operator.itemgetter(1))[:2]
            if most_common_allele[1]/float(depth_cnt) >= 0.45 and second_common_allele[1]/float(depth_cnt) >= 0.45:
                if alt != most_common_allele[0] and alt != second_common_allele[0]:
                    is_triple_alleles = True

    return (nlod, is_zero_coverage_in_normal, is_exist_in_normal, is_triple_alleles)


# -------------------------------------------------------------------------------------
# find median value of list.
# -------------------------------------------------------------------------------------
def median(list):
    quotient, remainder = divmod(len(list), 2)
    if remainder:
        return sorted(list)[quotient]
    return sum(sorted(list)[quotient - 1:quotient + 1]) / 2.


# -------------------------------------------------------------------------------------
# calculate active score of snp and indels.
# -------------------------------------------------------------------------------------
def cal_active_score(bamfile, region, refseq):
    samfile = pysam.AlignmentFile(bamfile, "rb")
    fasta = pysam.FastaFile(refseq)
    pos_to_score = defaultdict(tuple)
    chrom, start, end = region
    format_region = '_'.join(region)
    bq = 30
    # all sites in region
    for reads in samfile.pileup(chrom, int(start), int(end), truncate=True, stepper='nofilter'):
        bc_to_bq = defaultdict(list)
        bc_to_base = defaultdict(list)
        del_bc_to_bq = defaultdict(list)
        ins_bc_to_bq = defaultdict(list)
        bc_to_count = defaultdict(int)
        bc_to_del = defaultdict(int)
        bc_to_ins = defaultdict(int)
        ref = 0
        non_ref = 0
        delete = 0
        non_delete = 0
        insert = 0
        non_insert = 0
        if len(reads) < 100:
            continue
        pos_genome = reads.reference_pos + 1
        refbase = fasta.fetch(reference=format_region, start=pos_genome - int(start)-1,
                              end=pos_genome - int(start))
        refbase = refbase.upper()
        # all reads in a site
        for read in reads.pileups:
            qname = read.alignment.query_name
            bc = qname.split(':')[-2].split('-')[-1]
            #delete
            if read.query_position is None:
                del_bc_to_bq[bc].append(bq)
                bc_to_del[bc] += 1
                bc_to_count[bc] += 1
            #insert
            elif read.indel > 0:
                bq = read.alignment.query_qualities[read.query_position]
                ins_bc_to_bq[bc].append(bq)
                bc_to_ins[bc] += 1
                bc_to_count[bc] += 1
            else:
                if read.indel == 0:
                    base = read.alignment.query_sequence[read.query_position]
                    bq = read.alignment.query_qualities[read.query_position]
                    bc_to_base[bc].append(base)
                    bc_to_bq[bc].append(bq)
                    bc_to_count[bc] += 1
        # calculate score
        for bc, value in bc_to_ins.items():
            if value/bc_to_count[bc] > 0.9:
                bq = median(ins_bc_to_bq[bc])
                insert += bq
            else:
                non_insert += bq
        if non_insert == 0:
            non_insert = 2
        score_ins = insert / float(non_insert)

        for bc, value in bc_to_del.items():
            if value/bc_to_count[bc] > 0.9:
                bq = median(del_bc_to_bq[bc])
                delete += bq
            else:
                non_delete += bq
        if non_delete == 0:
            non_delete = 2
        score_del = delete / float(non_delete)

        for bc, base in bc_to_base.items():
            bq = median(bc_to_bq[bc])
            base = Counter(base).most_common(1)[0][0]
            if base == refbase:
                ref += bq
            else:
                non_ref += bq
        if ref == 0:
            ref = 2
        score_snp = non_ref / float(ref)
        raw_score = score_snp + score_del + score_ins
        pos_to_score[pos_genome] = (raw_score,score_del)

    return (region, pos_to_score)


# -------------------------------------------------------------------------------------
# convert variant type, reference base, variant base to output format
# -------------------------------------------------------------------------------------
def convert_to_vcf(orig_ref, orig_alt):
    vtype = '.'
    ref = orig_ref
    alt = orig_alt
    if len(orig_alt) == 1:
        vtype = 'SNP'
    elif orig_alt == 'DEL':
        vtype = 'SDEL'
    else:
        vals = orig_alt.split('|')
        if vals[0] in ('DEL', 'INS'):
            vtype = 'INDEL'
            ref = vals[1]
            alt = vals[2]

    return (ref, alt, vtype)


# -------------------------------------------------------------------------------------
# check if a locus is within or flanked by homopolymer region and/or low complexity region
# -------------------------------------------------------------------------------------
def is_hom_low_comp(chrom, pos, length, refb, altb, refgenome):
   # get reference bases for interval [pos-length, pos+length]
    refs = pysam.FastaFile(refgenome)
    chrom_length = refs.get_reference_length(chrom)
    pos0 = int(pos) - 1
    l_seq = refs.fetch(reference=chrom, start=max(0, pos0 - length), end=pos0).upper()
    r_seq_ref = refs.fetch(reference=chrom, start=pos0 + len(refb),
                          end=min(pos0 + len(refb) + length, chrom_length)).upper()
    r_seq_alt = refs.fetch(reference=chrom, start=pos0 + len(altb),
                          end=min(pos0 + len(altb) + length, chrom_length)).upper()
    ref_seq = l_seq + refb + r_seq_ref
    alt_seq = l_seq + altb + r_seq_alt
    # check homopolymer
    homo_a = ref_seq.find('A' * length) >= 0 or alt_seq.find('A' * length) >= 0
    homo_t = ref_seq.find('T' * length) >= 0 or alt_seq.find('T' * length) >= 0
    homo_g = ref_seq.find('G' * length) >= 0 or alt_seq.find('G' * length) >= 0
    homo_c = ref_seq.find('C' * length) >= 0 or alt_seq.find('C' * length) >= 0
    homo_p = homo_a or homo_t or homo_g or homo_c

    # check low complexity -- window length is 2 * homopolymer region. If any 2 nucleotide >= 99%
    len2 = 2 * length
    l_seq_lc = refs.fetch(reference=chrom, start=max(0, pos0 - len2), end=pos0).upper()
    r_seq_ref_lc = refs.fetch(reference=chrom, start=pos0 + len(refb),
                            end=min(pos0 + len(refb) + len2, chrom_length)).upper()  # ref seq
    r_seq_alt_lc = refs.fetch(reference=chrom, start=pos0 + len(altb),
                            end=min(pos0 + len(altb) + len2, chrom_length)).upper()  # alt seq
    ref_seq_lc = l_seq_lc + refb + r_seq_ref_lc
    alt_seq_lc = l_seq_lc + altb + r_seq_alt_lc
    low_comp = False

    # Ref seq
    total_len = len(ref_seq_lc)
    for i in range(total_len - len2):
        subseq = ref_seq_lc[i:(i + len2)]
        count_a = subseq.count('A')
        count_t = subseq.count('T')
        count_g = subseq.count('G')
        count_c = subseq.count('C')
        sorted_counts = sorted([count_a, count_t, count_g, count_c], reverse=True)
        top2_freq = 1.0 * (sorted_counts[0] + sorted_counts[1]) / len2
        if top2_freq >= 0.99:
            low_comp = True
            break

    # If ref seq is not LC, check alt seq
    if not low_comp:
        total_len = len(alt_seq_lc)
        for i in range(total_len - len2):
            subseq = alt_seq_lc[i:(i + len2)]
            count_a = subseq.count('A')
            count_t = subseq.count('T')
            count_g = subseq.count('G')
            count_c = subseq.count('C')
            sorted_counts = sorted([count_a, count_t, count_g, count_c], reverse=True)
            top2_freq = 1.0 * (sorted_counts[0] + sorted_counts[1]) / len2
            if top2_freq >= 0.99:
                low_comp = True
                break

    return (homo_p, low_comp)


# -------------------------------------------------------------------------------------
# filter variants
# -------------------------------------------------------------------------------------
def filter_variants(ref, alt, vtype, orig_alt, orig_ref, used_mt, strong_mt_cnt, chrom, pos, hp_len, ref_genome, mt_cnt,
                   allele_cnt, cvg, discord_pair_cnt, concord_pair_cnt, reverse_cnt, forward_cnt, lowq_reads, r1bc_end_pos,
                   r2bc_end_pos, r2primer_end_pos, primer_dist, dis_start_list, dis_end_list):
    # init output string
    fltr = ';'

    # low coverage filter
    if used_mt < 5:
        fltr += 'LM;'

    # check region for homopolymer or low complexity
    (is_homopolymer, is_low_complexity) = is_hom_low_comp(chrom, pos, hp_len, ref, alt, ref_genome)

    # homopolymer filter
    if is_homopolymer and 1.0 * mt_cnt[orig_alt] / used_mt < 0.99:
        fltr += 'HP;'

    # low complexity filter
    if is_low_complexity and 1.0 * mt_cnt[orig_alt] / used_mt < 0.99:
        fltr += 'LowC;'

    # strand bias and discordant pairs filter
    af_alt = 100.0 * allele_cnt[orig_alt] / cvg
    pairs = discord_pair_cnt[orig_alt] + concord_pair_cnt[orig_alt]  # total number of paired reads covering the pos
    if pairs >= 1000 and 1.0 * discord_pair_cnt[orig_alt] / pairs >= 0.5:
        fltr += 'DP;'
    elif af_alt <= 60.0:
        ref_r = reverse_cnt[orig_ref]
        ref_f = forward_cnt[orig_ref]
        alt_r = reverse_cnt[orig_alt]
        alt_f = forward_cnt[orig_alt]
        fisher = scipy.stats.fisher_exact([[ref_r, ref_f], [alt_r, alt_f]])
        odds_ratio = fisher[0]
        pvalue = fisher[1]
        if pvalue < 0.00001 and (odds_ratio >= 50 or odds_ratio <= 1.0 / 50):
            fltr += 'SB;'

    # base quality filter. Reject if more than 40% reads are lowQ
    if vtype == 'SNP' and orig_alt in allele_cnt.keys() and orig_alt in lowq_reads.keys():
        bq_alt = 1.0 * lowq_reads[orig_alt] / allele_cnt[orig_alt]
    else:
        bq_alt = 0.0
    if bq_alt > 0.4:
        fltr += 'LowQ;'

    # random end and fixed end position filters
    if vtype == 'SNP':
        # random end position filter
        end_base = 20  # distance to barcode end of the read
        # R1
        ref_le_end = sum(d <= end_base for d in r1bc_end_pos[orig_ref])  # number of REF R2 reads with distance <= endBase
        ref_gt_end = len(r1bc_end_pos[orig_ref]) - ref_le_end  # number of REF R2 reads with distance > endBase
        alt_le_end = sum(d <= end_base for d in r1bc_end_pos[orig_alt])  # number of ALT R2 reads with distance <= endBase
        alt_gt_end = len(r1bc_end_pos[orig_alt]) - alt_le_end  # number of ALT R2 reads with distance > endBase
        fisher = scipy.stats.fisher_exact([[ref_le_end, ref_gt_end], [alt_le_end, alt_gt_end]])
        odds_ratio = fisher[0]
        pvalue = fisher[1]
        if pvalue < 0.001 and odds_ratio < 0.05 and af_alt <= 60.0:
            fltr += 'R1CP;'
        # R2
        ref_le_end = sum(d <= end_base for d in r2bc_end_pos[orig_ref])  # number of REF R2 reads with distance <= endBase
        ref_gt_end = len(r2bc_end_pos[orig_ref]) - ref_le_end  # number of REF R2 reads with distance > endBase
        alt_le_end = sum(d <= end_base for d in r2bc_end_pos[orig_alt])  # number of ALT R2 reads with distance <= endBase
        alt_gt_end = len(r2bc_end_pos[orig_alt]) - alt_le_end  # number of ALT R2 reads with distance > endBase
        fisher = scipy.stats.fisher_exact([[ref_le_end, ref_gt_end], [alt_le_end, alt_gt_end]])
        odds_ratio = fisher[0]
        pvalue = fisher[1]
        if pvalue < 0.001 and odds_ratio < 0.05 and af_alt <= 60.0:
            fltr += 'R2CP;'

        # fixed end position filter
        end_base = primer_dist  # distance to primer end of the read
        ref_le_end = sum(d <= end_base for d in r2primer_end_pos[orig_ref])  # number of REF R2 reads with distance <= endBase
        ref_gt_end = len(r2primer_end_pos[orig_ref]) - ref_le_end  # number of REF R2 reads with distance > endBase
        alt_le_end = sum(d <= end_base for d in r2primer_end_pos[orig_alt])  # number of ALT R2 reads with distance <= endBase
        alt_gt_end = len(r2primer_end_pos[orig_alt]) - alt_le_end  # number of ALT R2 reads with distance > endBase
        fisher = scipy.stats.fisher_exact([[ref_le_end, ref_gt_end], [alt_le_end, alt_gt_end]])
        odds_ratio = fisher[0]
        pvalue = fisher[1]
        # reject if variant is clustered within 2 bases from primer sequence due to possible enzyme initiation error
        if alt_le_end + alt_gt_end > 0:
            if 1.0 * alt_le_end / (alt_le_end + alt_gt_end) >= 0.98 or (pvalue < 0.001 and odds_ratio < 1.0 / 20):
                fltr += 'PrimerCP;'

    median_start = median(dis_start_list)
    median_end = median(dis_end_list)
    if median_start <= 10:
        deviation_start = [abs(start-median_start) for start in dis_start_list]
        if median(deviation_start) <= 3:
            fltr += 'CPos'
    if median_end <= 10:
        deviation_end = [abs(end-median_end) for end in dis_end_list]
        if median(deviation_end) <= 3:
            fltr += 'CPos'

    # done
    return fltr


# -------------------------------------------------------------------------------------
# function to call variants
# -------------------------------------------------------------------------------------
def vc(bamfile, chrom, pos, min_bq, min_mq, rpb, hp_len, mismatch_thr, mt_drop, max_mt, tlod_threshold, primer_dist, refgenome):
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    idx = 0
    cvg = 0
    bc_dict = defaultdict(lambda: defaultdict(list))
    allbc_dict = defaultdict(list)
    allele_cnt = defaultdict(int)
    mt_cnt = defaultdict(int)
    r1_bc_end_pos = defaultdict(list)
    r2_bc_end_pos = defaultdict(list)
    r2_primer_end_pos = defaultdict(list)
    mt3_cnt = 0
    mt5_cnt = 0
    mt7_cnt = 0
    mt10_cnt = 0
    strong_mt_cnt = defaultdict(int)
    pred_index = defaultdict(lambda: defaultdict(float))
    final_dict = defaultdict(float)
    r1_cnt = defaultdict(int)
    r2_cnt = defaultdict(int)
    forward_cnt = defaultdict(int)
    reverse_cnt = defaultdict(int)
    concord_pair_cnt = defaultdict(int)
    discord_pair_cnt = defaultdict(int)
    mismatch_cnt = defaultdict(float)
    bq_sum = defaultdict(int)
    low_q_reads = defaultdict(int)
    allele_freq = defaultdict(float)
    likelihood_m0 = 1.0
    likelihood_somatic = 1.0
    divide = 1.0
    nearby_indel_cnt = defaultdict(int)
    dis_start_list = []
    dis_end_list = []

    # set threshold for strongMT based on mtDepth
    if rpb < 1.5:
        smt = 2.0
    elif rpb < 3.0:
        smt = 3.0
    else:
        smt = 4.0

    # get reference base
    refseq = pysam.FastaFile(refgenome)
    orig_ref = refseq.fetch(reference=chrom, start=int(pos) - 1, end=int(pos))
    orig_ref = orig_ref.upper()

    # pile up reads
    for read in samfile.pileup(region=chrom + ':' + pos + ':' + pos, truncate=True, max_depth=1000000,
                              stepper='nofilter'):
        for pileup_read in read.pileups:
            # read ID
            qname = pileup_read.alignment.query_name
            qname_split = qname.split(":")
            readid = ':'.join(qname_split[:-2])
            # barcode sequence
            bc = qname_split[-2]
            # duplex tag - temporary hack from end of readid - should be CC, TT, or NN for duplex runs
            duplexTag = qname_split[-3]
            # mapping quality
            mq = pileup_read.alignment.mapping_quality
            # get NM tag
            nm = 0
            all_tags = pileup_read.alignment.tags
            for (tag, value) in all_tags:
                if tag == 'NM':
                    nm = value
                    break
            # count number of INDELs in the read sequence
            n_indel = 0
            cigar = pileup_read.alignment.cigar
            cigar_order = 1
            left_sp = 0  # soft clipped bases on the left
            right_sp = 0  # soft clipped bases on the right
            for (op, value) in cigar:
                # 1 for insertion
                if op == 1 or op == 2:
                    n_indel += value
                if cigar_order == 1 and op == 4:
                    left_sp = value
                if cigar_order > 1 and op == 4:
                    right_sp += value
                cigar_order += 1

            # Number of mismatches except INDEL, including softcilpped sequences
            mismatch = max(0, nm - n_indel)
            # read length, including softclip
            read_len = pileup_read.alignment.query_length
            # calculate mismatch per 100 bases
            mismatch_per_100b = 100.0 * mismatch / read_len if read_len > 0 else 0.0

            # paired read
            if pileup_read.alignment.is_read1:
                pair_order = 'R1'
            if pileup_read.alignment.is_read2:
                pair_order = 'R2'

            # +/- strand
            strand = 'Reverse' if pileup_read.alignment.is_reverse else 'Forward'

            # coverage -- read, not fragment
            cvg += 1

            if cigar[0][0] == 4:
                true_start = pileup_read.alignment.reference_start - cigar[0][1]
            else:
                true_start = pileup_read.alignment.reference_start
            if cigar[-1][0] == 4:
                true_end = pileup_read.alignment.reference_end + cigar[-1][1]
            else:
                true_end = pileup_read.alignment.reference_end
            if pileup_read.alignment.is_reverse:
                dis_start = true_end - int(pos)
                dis_end = int(pos) - true_start
            else:
                dis_start = int(pos) - true_start
                dis_end = true_end - int(pos)
            dis_start_list.append(dis_start)
            dis_end_list.append(dis_end)

            # check if the site is the beginning of insertion
            if pileup_read.indel > 0:
                site = pileup_read.alignment.query_sequence[pileup_read.query_position]
                inserted = pileup_read.alignment.query_sequence[
                            (pileup_read.query_position + 1): (pileup_read.query_position + 1 + pileup_read.indel)]
                base = 'INS|' + site + '|' + site + inserted
                bq = pileup_read.alignment.query_qualities[pileup_read.query_position]
                bq_sum[base] += bq
                # inclusion condition
                inc_cond = bq >= min_bq and mq >= min_mq and mismatch_per_100b <= mismatch_thr
                allele_cnt[base] += 1
                mismatch_cnt[base] += mismatch_per_100b
                if pair_order == 'R1':
                    r1_cnt[base] += 1
                if pair_order == 'R2':
                    r2_cnt[base] += 1

                if strand == 'Reverse':
                    reverse_cnt[base] += 1
                else:
                    forward_cnt[base] += 1

                for p in range(int(pos)-5, int(pos)):
                    nearby_indel_cnt[p] += 1
                for p in range(int(pos)+1, int(pos)+6):
                    nearby_indel_cnt[p] += 1

            # check if the site is the beginning of deletion
            elif pileup_read.indel < 0:
                site = pileup_read.alignment.query_sequence[pileup_read.query_position]
                deleted = refseq.fetch(reference=chrom, start=int(pos), end=int(pos) + abs(pileup_read.indel))
                deleted = deleted.upper()
                base = 'DEL|' + site + deleted + '|' + site
                bq = pileup_read.alignment.query_qualities[pileup_read.query_position]
                bq_sum[base] += bq
                # inclusion condition
                inc_cond = bq >= min_bq and mq >= min_mq and mismatch_per_100b <= mismatch_thr
                allele_cnt[base] += 1
                mismatch_cnt[base] += mismatch_per_100b
                if pair_order == 'R1':
                    r1_cnt[base] += 1
                if pair_order == 'R2':
                    r2_cnt[base] += 1

                if strand == 'Reverse':
                    reverse_cnt[base] += 1
                else:
                    forward_cnt[base] += 1

                dist_del = abs(pileup_read.indel)
                for p in range(int(pos)-5, int(pos)):
                    nearby_indel_cnt[p] += 1
                for p in range(int(pos)+1+dist_del, int(pos)+6+dist_del):
                    nearby_indel_cnt[p] += 1

            # site is not beginning of any INDEL
            else:
                # If the site ifself is a deletion, set quality = minBQ
                if pileup_read.is_del:
                    base = 'DEL'
                    bq = min_bq
                    bq_sum[base] += bq
                    # inclusion condition
                    inc_cond = bq >= min_bq and mq >= min_mq and mismatch_per_100b <= mismatch_thr
                # if the site is a regular locus,
                else:
                    base = pileup_read.alignment.query_sequence[
                        pileup_read.query_position]  # note: query_sequence includes soft clipped bases
                    bq = pileup_read.alignment.query_qualities[pileup_read.query_position]
                    bq_sum[base] += bq
                    # count the number of low quality reads (less than Q20 by default) for each base
                    if bq < min_bq:
                        low_q_reads[base] += 1
                    # inclusion condition
                    inc_cond = bq >= min_bq and mq >= min_mq and mismatch_per_100b <= mismatch_thr
                    if pair_order == 'R1':
                        # distance to the barcode end in R1;
                        if pileup_read.alignment.is_reverse:
                            dist_to_bc_end = pileup_read.alignment.query_alignment_length - (pileup_read.query_position - left_sp)
                        else:
                            dist_to_bc_end = pileup_read.query_position - left_sp
                        if inc_cond:
                            r1_bc_end_pos[base].append(dist_to_bc_end)
                            r1_cnt[base] += 1
                    if pair_order == 'R2':
                        # distance to the barcode and/or primer end in R2. Different cases for forward and reverse strand
                        if pileup_read.alignment.is_reverse:
                            dist_to_bc_end = pileup_read.query_position - left_sp
                            dist_to_primer_end = pileup_read.alignment.query_alignment_length - (
                                        pileup_read.query_position - left_sp)
                        else:
                            dist_to_bc_end = pileup_read.alignment.query_alignment_length - (pileup_read.query_position - left_sp)
                            dist_to_primer_end = pileup_read.query_position - left_sp
                        if inc_cond:
                            r2_bc_end_pos[base].append(dist_to_bc_end)
                            r2_primer_end_pos[base].append(dist_to_primer_end)
                        r2_cnt[base] += 1

                    if strand == 'Reverse':
                        reverse_cnt[base] += 1
                    else:
                        forward_cnt[base] += 1

                allele_cnt[base] += 1
                mismatch_cnt[base] += mismatch_per_100b

            # count total number of fragments and MTs
            if readid not in allbc_dict[bc]:
                allbc_dict[bc].append(readid)

            # decide which read goes into analysis
            if inc_cond:
                if readid not in bc_dict[bc]:
                    prob = pow(10.0, -bq / 10.0)
                    readinfo = [base, prob, pair_order]
                    bc_dict[bc][readid].append(readinfo)
                elif base == bc_dict[bc][readid][0][0] or base in ['N', '*']:
                    bc_dict[bc][readid][0][1] = max((pow(10.0, -bq / 10.0), bc_dict[bc][readid][0][1]))
                    bc_dict[bc][readid][0][2] = 'Paired'
                    if base == bc_dict[bc][readid][0][0]:
                        concord_pair_cnt[base] += 1
                else:
                    del bc_dict[bc][readid]
                    discord_pair_cnt[base] += 1
                #drop N, after cluster UMI.
                if base in ['N', '*']:
                    del bc_dict[bc]

    # total number of MT, fragments, reads, including those dropped from analysis
    all_mt = len(allbc_dict)
    all_frag = sum([len(allbc_dict[bc]) for bc in allbc_dict])


    # MTs used
    used_mt = len(bc_dict)

    # done if zero coverage (note hack for 41 blank output fields!)
    if used_mt == 0:
        out_long = '\t'.join([chrom, pos, orig_ref] + [''] * 41 + ['Zero_Coverage'])
        return out_long

    bc_keys = bc_dict.keys()
    used_frag = sum([len(bc_dict[bc]) for bc in bc_keys])

    total_r1 = sum(r1_cnt.values())
    total_r2 = sum(r2_cnt.values())

    for allele, cnt in allele_cnt.items():
        allele_freq[allele] = float(cnt)/float(cvg)
    if len(allele_freq) != 1:
        sorted_list = sorted(allele_freq.items(), key=operator.itemgetter(1), reverse=True)
        max_base, max_freq = sorted_list[0]
        second_base, second_freq = sorted_list[1]
        orig_alt = second_base if max_base == orig_ref else max_base
        alt_freq = second_freq if max_base == orig_ref else max_freq

        # calculate tlod score and determined variants
        for bc in bc_keys:
            p0_right_bc, p1_right_bc, base = cal_prob(bc_dict[bc], orig_ref, orig_alt, alt_freq)
            likelihood_m0 *= p0_right_bc
            likelihood_somatic *= p1_right_bc
            divide *= p1_right_bc/p0_right_bc
            mt_cnt[base] += 1

            if len(bc_dict[bc]) >= 3:
                mt3_cnt += 1
            if len(bc_dict[bc]) >= 5:
                mt5_cnt += 1
            if len(bc_dict[bc]) >= 7:
                mt7_cnt += 1
            if len(bc_dict[bc]) >= 10:
                mt10_cnt += 1

        try:
            tlod = math.log10(divide)
        except ValueError:
            tlod = 0

        if tlod >= tlod_threshold:
            alt_pi = 5
            second_pi = 5
            # convert from internal smCounter format to format needed for output
            (ref, alt, vtype) = convert_to_vcf(orig_ref, orig_alt)

            # apply filters if PI >= 5 (at least 2 MTs), and locus not in a deletion
            fltr = ';'
            if alt_pi >= 5 and vtype in ('SNP', 'INDEL'):
                fltr = filter_variants(ref, alt, vtype, orig_alt, orig_ref, used_mt, strong_mt_cnt, chrom, pos, hp_len, refgenome, mt_cnt,
                                      allele_cnt, cvg, discord_pair_cnt, concord_pair_cnt, reverse_cnt, forward_cnt, low_q_reads,
                                      r1_bc_end_pos, r2_bc_end_pos, r2_primer_end_pos, primer_dist, dis_start_list, dis_end_list)

            # identify possible bi-allelic variants - top 2 alleles are non-reference and both VMFs >= 45%. Not necessarily passing the filters.
            mf_alt = 1.0 * mt_cnt[max_base] / used_mt  # MT fraction of the base with the highest PI
            mf_alt2 = 1.0 * mt_cnt[second_base] / used_mt  # MT fraction of the base with the second highest PI
            if max_base != orig_ref and second_base != orig_ref and mf_alt >= 0.45 and mf_alt2 >= 0.45:  # conditions to be considered bi-allelic

                # convert from internal smCounter format to format needed for output
                orig_alt2 = second_base
                (ref2, alt2, vtype2) = convert_to_vcf(orig_ref, orig_alt)

                # apply filters to 2nd variant if PI2 >= 5 (at least 2 MTs), and locus not in a deletion
                fltr2 = ';'
                if second_pi >= 5 and vtype2 in ('SNP', 'INDEL'):
                    fltr2 = filter_variants(ref2, alt2, vtype2, orig_alt2, orig_ref, used_mt, strong_mt_cnt, chrom, pos, hp_len,
                                           refgenome, mt_cnt, allele_cnt, cvg, discord_pair_cnt, concord_pair_cnt, reverse_cnt,
                                           forward_cnt, low_q_reads, r1_bc_end_pos, r2_bc_end_pos, r2_primer_end_pos, primer_dist)

                # prepare output for bi-allelic variants (if var2 is filtered, regardless of var1, do nothing. output var1 only)
                if fltr == ';' and fltr2 == ';':  # both var1 and var2 pass the filters -- this is a bi-allelic variant. var1's statistics (MT, DP, etc) are reported
                    alt = alt + ',' + alt2
                    vtype = vtype.lower() + ',' + vtype2.lower()
                elif fltr != ';' and fltr2 == ';':  # if var1 is filtered and the var2 passes, then it's a single variant of var2
                    alt = alt2
                    fltr = fltr2
                    orig_alt = orig_alt2

            # build detailed output vector
            frac_alt = round((1.0 * allele_cnt[orig_alt] / cvg), 4)  # based on all reads, including the excluded reads
            frac_a = round((1.0 * allele_cnt['A'] / cvg), 4)
            frac_t = round((1.0 * allele_cnt['T'] / cvg), 4)
            frac_g = round((1.0 * allele_cnt['G'] / cvg), 4)
            frac_c = round((1.0 * allele_cnt['C'] / cvg), 4)
            fracs = (allele_cnt['A'], allele_cnt['T'], allele_cnt['G'], allele_cnt['C'], frac_a, frac_t, frac_g, frac_c)

            mt_f_alt = round((1.0 * mt_cnt[orig_alt] / used_mt), 4)  # based on only used MTs
            mt_f_a = round((1.0 * mt_cnt['A'] / used_mt), 4)
            mt_f_t = round((1.0 * mt_cnt['T'] / used_mt), 4)
            mt_f_g = round((1.0 * mt_cnt['G'] / used_mt), 4)
            mt_f_c = round((1.0 * mt_cnt['C'] / used_mt), 4)
            mts = (
            mt3_cnt, mt5_cnt, mt7_cnt, mt10_cnt, mt_cnt['A'], mt_cnt['T'], mt_cnt['G'], mt_cnt['C'], mt_f_a, mt_f_t, mt_f_g, mt_f_c)

            strong_mt = (strong_mt_cnt['A'], strong_mt_cnt['T'], strong_mt_cnt['G'], strong_mt_cnt['C'])
            pred_idx = (round(final_dict['A'], 2), round(final_dict['T'], 2), round(final_dict['G'], 2), round(final_dict['C'], 2))

            outvec = [chrom, pos, ref, alt, orig_alt, vtype, cvg, all_frag, all_mt, used_frag, used_mt, round(final_dict[orig_alt], 2),
                     allele_cnt[orig_alt], frac_alt, mt_cnt[orig_alt], mt_f_alt, strong_mt_cnt[orig_alt], round(tlod, 2)]
            outvec.extend(fracs)
            outvec.extend(mts)
            outvec.extend(strong_mt)
            outvec.extend(pred_idx)
            outvec.append(fltr)
            out_long = '\t'.join((str(x) for x in outvec))

            return (out_long, nearby_indel_cnt)


# ------------------------------------------------------------------------------------------------
# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace
# ------------------------------------------------------------------------------------------------
def vc_wrapper(*args):
    output = vc(*args)

    return output


# ------------------------------------------------------------------------------------------------
# global for argument parsing (hack that works when calling from either command line or pipeline)
# ------------------------------------------------------------------------------------------------
parser = None


def argParseInit():  # this is done inside a function because multiprocessing module imports the script
    global parser
    parser = argparse.ArgumentParser(description='Variant calling using molecular barcodes', fromfile_prefix_chars='@')
    parser.add_argument('--outPrefix', default=None, required=True, help='prefix for output files')
    parser.add_argument('--bamFile', default=None, required=True, help='BAM file')
    parser.add_argument('--normal_bamfile', default=None, required=True, help='normal sample BAM file')
    parser.add_argument('--bedTarget', default=None, required=True, help='BED file for target region')
    parser.add_argument('--fastaTarget', default=None, required=True, help='FASTA file for target region')
    parser.add_argument('--rpb', default=None, required=True, type=float, help='Mean read pairs per MT')
    parser.add_argument('--nCPU', type=int, default=1, help='number of CPUs to use in parallel')
    parser.add_argument('--minBQ', type=int, default=20, help='minimum base quality allowed for analysis')
    parser.add_argument('--minMQ', type=int, default=30, help='minimum mapping quality allowed for analysis')
    parser.add_argument('--hpLen', type=int, default=10, help='Minimum length for homopolymers')
    parser.add_argument('--mismatchThr', type=float, default=6.0,
                       help='average number of mismatches per 100 bases allowed')
    parser.add_argument('--mtDrop', type=int, default=0, help='Drop MTs with lower than or equal to X reads.')
    parser.add_argument('--maxMT', type=int, default=0,
                       help='Randomly downsample to X MTs (max number of MTs at any position). If set to 0 (default), maxMT = 2.0 * mean MT depth')
    parser.add_argument('--primerDist', type=int, default=2, help='filter variants that are within X bases to primer')
    parser.add_argument('--threshold', type=int, default=0,
                       help='Minimum prediction index for a variant to be called. Must be non-negative. Typically ranges from 10 to 60. If set to 0 (default), smCounter will choose the appropriate cutoff based on the mean MT depth.')
    parser.add_argument('--min_frequency', type=float, default=0.01, help='Minimum frequency for variants')
    parser.add_argument('--min_active_score', type=float, default=0.01, help='Minimum active score for positions')
    parser.add_argument('--tlod_threshold', type=float, default=6.3, help='Minimum tlod for variants')
    parser.add_argument('--nlod_threshold', type=float, default=10, help='Minimum nlod for variants')
    parser.add_argument('--refGenome', default='/qgen/home/rvijaya/downloads/alt_hap_masked_ref/ucsc.hg19.fasta')
    parser.add_argument('--bedTandemRepeats', default='/qgen/home/xuc/UCSC/simpleRepeat.bed',
                       help='bed for UCSC tandem repeats')
    parser.add_argument('--bedRepeatMaskerSubset', default='/qgen/home/xuc/UCSC/SR_LC_SL.nochr.bed',
                       help='bed for RepeatMasker simple repeats, low complexity, microsatellite regions')
    parser.add_argument('--bedtoolsPath', default='/qgen/bin/bedtools-2.25.0/bin/', help='path to bedtools')
    parser.add_argument('--runPath', default=None, help='path to working directory')
    parser.add_argument('--logFile', default=None, help='log file')
    parser.add_argument('--paramFile', default=None,
                       help='optional parameter file that contains the above paramters. if specified, this must be the only parameter, except for --logFile.')


# --------------------------------------------------------------------------------------
# main function
# --------------------------------------------------------------------------------------
def main(args):
    # log run start
    time_start = datetime.datetime.now()
    print("smCounter started at " + str(time_start))

    # if argument parser global not assigned yet, initialize it
    if parser == None:
      argParseInit()

    # get arguments passed in via a lambda object (e.g. from upstream pipeline)
    if type(args) is not argparse.Namespace:
        args_list = []
        for arg_name, arg_val in args.iteritems():
            args_list.append("--{0}={1}".format(arg_name, arg_val))
        args = parser.parse_args(args_list)

    # get arguments from disk file specified on command line (warning: this silently deletes all actual command line parameters)
    elif args.paramFile != None:
        args = parser.parse_args(("@" + args.paramFile,))

    # echo all parameters to the log file
    for arg_name, arg_val in vars(args).items():
        print(arg_name, arg_val)

    # change working directory to runDir
    if args.runPath != None:
        os.chdir(args.runPath)

    # make list of loci to call variants
    loc_list = []
    region_to_pos = defaultdict(dict)
    bamfile = args.bamFile
    refseq = args.fastaTarget
    raw_regions = defaultdict(str)
    for seq in open(refseq, 'r'):
        if seq.startswith('>'):
            pos = tuple(seq.strip()[1:].split('_'))
            raw_regions[pos] = ''
        else:
            raw_regions[pos] += seq.strip()
    pool1 = multiprocessing.Pool(processes=args.nCPU)
    result1 = [pool1.apply_async(cal_active_score, args=(bamfile, region, refseq)) for region in raw_regions]
    for p1 in result1:
        region_to_pos[p1.get()[0]] = p1.get()[1]
    pool1.close()
    pool1.join()

    # for region in raw_regions:
    #     nnn = cal_active_score(bamfile, region, refseq)
    #     region_to_pos[nnn[0]] = nnn[1]

    for region,pos_to_score in region_to_pos.items():
        pos_to_score = sorted(pos_to_score.items())
        for pos,score in pos_to_score:
            if score[0] >= args.min_active_score:
                loc_list.append((region[0],str(pos)))
                if score[1] >= args.min_active_score:
                    loc_list.append((region[0],str(pos-1)))
                    loc_list.append((region[0],str(pos+1)))
    loc_list = sorted(list(set(loc_list)))

    # call variants in parallel
    pool = multiprocessing.Pool(processes=args.nCPU)
    results = [pool.apply_async(vc_wrapper, args=(
    args.bamFile, x[0], x[1], args.minBQ, args.minMQ, args.rpb, args.hpLen, args.mismatchThr, args.mtDrop,
    args.maxMT, args.tlod_threshold, args.primerDist, args.refGenome)) for x in loc_list]
    options = [p.get() for p in results]
    pool.close()
    pool.join()
    # options = []
    output = []
    # for x in loc_list:
    #     options.append(vc_wrapper(args.bamFile, x[0], x[1], args.minBQ, args.minMQ, args.mtDepth, args.rpb, args.hpLen, args.mismatchThr, args.mtDrop,
    # args.maxMT, args.primerDist, args.refGenome))
    options = [line for line in options if line!= None]
    nearby_indel_cnt = defaultdict(int)
    for out in options:
        line, indel_cnt = out
        output.append(line)
        for pos, cnt in indel_cnt.items():
            nearby_indel_cnt[pos] += cnt

    # check for exceptions thrown by vc()
    for idx in range(len(output)):
        line = output[idx]
        if line.startswith("Exception thrown!"):
            print(line)
            raise Exception("Exception thrown in vc() at location: " + str(loc_list[idx]))

    # report start of variant filtering
    print("begin variant filtering and output")

    # merge and sort RepeatMasker tracks (could be done prior to run)  Note: assuming TRF repeat already merged and sorted!!
    bed_exe = args.bedtoolsPath
    bed_repeat_masker = args.outPrefix + '.tmp.repeatMasker.bed'
    subprocess.check_call(
        bed_exe + ' merge -c 4 -o distinct -i ' + args.bedRepeatMaskerSubset + ' | ' + bed_exe + ' sort -i - > ' + bed_repeat_masker,
        shell=True)

    # merge and sort target region
    bed_target = args.outPrefix + '.tmp.target.bed'
    subprocess.check_call(bed_exe + ' merge -i ' + args.bedTarget + ' | ' + bed_exe + ' sort -i - > ' + bed_target,
                          shell=True)

    # intersect 2 repeats tracks with target region
    subprocess.check_call(
        bed_exe + ' intersect -a ' + args.bedTandemRepeats + ' -b ' + bed_target + ' | ' + bed_exe + ' sort -i - > ' + args.outPrefix + '.tmp.target.repeats1.bed',
        shell=True)
    subprocess.check_call(
        bed_exe + ' intersect -a ' + bed_repeat_masker + ' -b ' + bed_target + ' | ' + bed_exe + ' sort -i - > ' + args.outPrefix + '.tmp.target.repeats2.bed',
        shell=True)

    # read in tandem repeat list
    trf_regions = defaultdict(list)
    for line in open(args.outPrefix + '.tmp.target.repeats1.bed', 'r'):
        vals = line.strip().split()
        (chrom, region_start, region_end) = vals[0:3]
        trf_regions[chrom].append((int(region_start), int(region_end), "RepT;"))

    # read in simple repeat, low complexity, satelite list
    rm_regions = defaultdict(list)
    for line in open(args.outPrefix + '.tmp.target.repeats2.bed', 'r'):
        (chrom, region_start, region_end, type_codes) = line.strip().split()
        rep_types = []
        for type_code in type_codes.split(","):
            if type_code == 'Simple_repeat':
                rep_types.append('RepS')
            elif type_code == 'Low_complexity':
                rep_types.append('LowC')
            elif type_code == 'Satellite':
                rep_types.append('SL')
            else:
                rep_types.append('Other_Repeat')
        rep_type = ";".join(rep_types) + ";"
        rm_regions[chrom].append((int(region_start), int(region_end), rep_types))

    # remove intermediate files
    os.remove(args.outPrefix + '.tmp.target.bed')
    os.remove(args.outPrefix + '.tmp.repeatMasker.bed')
    os.remove(args.outPrefix + '.tmp.target.repeats1.bed')
    os.remove(args.outPrefix + '.tmp.target.repeats2.bed')

    # set up header columns (Note: "headerAll" must parallel the output of the vc() function.)
    header_all = (
    'CHROM', 'POS', 'REF', 'ALT', 'ORIG_ALT', 'TYPE', 'DP', 'FR', 'MT', 'UFR', 'UMT', 'PI', 'VDP', 'VAF', 'VMT', 'VMF',
    'VSM', 'TLOD', 'DP_A', 'DP_T', 'DP_G', 'DP_C', 'AF_A', 'AF_T', 'AF_G', 'AF_C', 'MT_3RPM', 'MT_5RPM', 'MT_7RPM',
    'MT_10RPM', 'UMT_A', 'UMT_T', 'UMT_G', 'UMT_C', 'UMF_A', 'UMF_T', 'UMF_G', 'UMF_C', 'VSM_A', 'VSM_T', 'VSM_G',
    'VSM_C', 'PI_A', 'PI_T', 'PI_G', 'PI_C', 'FILTER')
    header_variants = (
    'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'DP', 'MT', 'UMT', 'PI', 'THR', 'VMT', 'VMF', 'VSM', 'TLOD', 'FILTER')

    # set up hash of variable fields
    header_all_index = {}
    for i in range(len(header_all)):
      header_all_index[header_all[i]] = i

    # ALL repeats filter. If MT fraction < 40% and the variant is inside the tandem repeat region, reject.
    for i in range(len(output)):
        outline = output[i]
        line_list = outline.split('\t')
        chrom_tr = line_list[header_all_index['CHROM']]
        alt_tr = line_list[header_all_index['ALT']]
        try:
            pos_tr = int(line_list[header_all_index['POS']])
        except ValueError:
            continue
        try:
            alt_mt_frac_tr = float(line_list[header_all_index['VMF']])
        except ValueError:
            continue
        try:
            pred = int(float(line_list[header_all_index['PI']]))
        except ValueError:
            pred = 0

        if pred >= 5 and alt_tr != 'DEL':
            # check tandem repeat from TRF if MT fraction < 40%
            if alt_mt_frac_tr < 40:
                for (loc_l, loc_r, rep_type) in trf_regions[chrom_tr]:
                    if loc_l < pos_tr <= loc_r:
                        line_list[-1] += rep_type
                        break

            # check simple repeat, lc, sl from RepeatMasker
            for (loc_l, loc_r, rep_type) in rm_regions[chrom_tr]:
                if loc_l < pos_tr <= loc_r:
                    line_list[-1] += rep_type
                    break

        line_list[-1] = 'PASS' if line_list[-1] == ';' else line_list[-1].strip(';')
        output[i] = '\t'.join(line_list)

    # VCF header
    header_vcf = \
        '##fileformat=VCFv4.2\n' + \
        '##reference=GRCh37\n' + \
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">\n' + \
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">\n' + \
        '##INFO=<ID=MT,Number=1,Type=Integer,Description="Total MT depth">\n' + \
        '##INFO=<ID=UMT,Number=1,Type=Integer,Description="Filtered MT depth">\n' + \
        '##INFO=<ID=PI,Number=1,Type=Float,Description="Variant prediction index">\n' + \
        '##INFO=<ID=THR,Number=1,Type=Integer,Description="Variant prediction index minimum threshold">\n' + \
        '##INFO=<ID=VMT,Number=1,Type=Integer,Description="Variant MT depth">\n' + \
        '##INFO=<ID=VMF,Number=1,Type=Float,Description="Variant MT fraction">\n' + \
        '##INFO=<ID=VSM,Number=1,Type=Integer,Description="Variant strong MT depth">\n' + \
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + \
        '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">\n' + \
        '##FORMAT=<ID=VF,Number=1,Type=Float,Description="Variant MT fraction, same as VMF">\n' + \
        '##FILTER=<ID=RepT,Description="Variant in simple tandem repeat region, as defined by Tandem Repeats Finder">\n' + \
        '##FILTER=<ID=RepS,Description="Variant in simple repeat region, as defined by RepeatMasker">\n' + \
        '##FILTER=<ID=LowC,Description="Variant in low complexity region, as defined by RepeatMasker">\n' + \
        '##FILTER=<ID=SL,Description="Variant in micro-satelite region, as defined by RepeatMasker">\n' + \
        '##FILTER=<ID=HP,Description="Inside or flanked by homopolymer region">\n' + \
        '##FILTER=<ID=LM,Description="Low coverage (fewer than 5 MTs)">\n' + \
        '##FILTER=<ID=LSM,Description="Fewer than 2 strong MTs">\n' + \
        '##FILTER=<ID=SB,Description="Strand bias">\n' + \
        '##FILTER=<ID=LowQ,Description="Low base quality (mean < 22)">\n' + \
        '##FILTER=<ID=MM,Description="Too many genome reference mismatches in reads (default threshold is 6.5 per 100 bases)">\n' + \
        '##FILTER=<ID=DP,Description="Too many discordant read pairs">\n' + \
        '##FILTER=<ID=R1CP,Description="Variants are clustered at the end of R1 reads">\n' + \
        '##FILTER=<ID=R2CP,Description="Variants are clustered at the end of R2 reads">\n' + \
        '##FILTER=<ID=PrimerCP,Description="Variants are clustered immediately after the primer, possible enzyme initiation error">\n' + \
        '\t'.join(('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.outPrefix)) + '\n'

    # set cutoff value for about 20 FP/Mb
    threshold = args.threshold

    # open output files
    out_all = open(args.outPrefix + '.smCounter.all.txt', 'w')
    out_variants = open(args.outPrefix + '.smCounter.somatic.cut.txt', 'w')
    out_vcf = open(args.outPrefix + '.smCounter.somatic.cut.vcf', 'w')

    # write column headers
    out_all.write('\t'.join(header_all) + '\n')
    out_variants.write('\t'.join(header_variants) + '\n')
    out_vcf.write(header_vcf)

    for line in output:
        # write to the detailed output
        out_all.write(line)
        out_all.write("\n")

        # unpack text fields
        fields = line.split('\t')

        # skip if no PI
        pi = fields[header_all_index['PI']]
        if len(pi) == 0:
            continue

        # get ALT and prediction index
        alt = fields[header_all_index['ALT']]
        orig_alt = fields[header_all_index['ORIG_ALT']]
        qual = str(int(float(pi)))  # truncate PI to conform to VCF phred-like tradition

        # write to vcf file and short output
        if int(qual) >= threshold and alt != 'DEL':  # if PI > threshold, write to vcf (regardless of filters)

            # parse fields needed from main data vector
            chrom = fields[header_all_index['CHROM']]
            pos = fields[header_all_index['POS']]
            ref = fields[header_all_index['REF']]
            vtype = fields[header_all_index['TYPE']]
            dp = fields[header_all_index['DP']]
            mt = fields[header_all_index['MT']]
            umt = fields[header_all_index['UMT']]
            vmt = fields[header_all_index['VMT']]
            vmf = fields[header_all_index['VMF']]
            vsm = fields[header_all_index['VSM']]
            tlod = fields[header_all_index['TLOD']]
            filter = fields[header_all_index['FILTER']]
            thr = str(threshold)
            info = ';'.join(('TYPE=' + vtype, 'DP=' + dp, 'MT=' + mt, 'UMT=' + umt, 'PI=' + pi, 'THR=' + thr, 'VMT=' + vmt,
                             'VMF=' + vmf, 'VSM=' + vsm, 'TLOD=' + tlod))

            if float(vmf) >= args.min_frequency:
                nlod, is_zero_coverage_in_normal, is_exist_in_normal, is_triple_alleles \
                    = variants_classification(chrom, pos, ref, orig_alt, args.normal_bamfile, args.minBQ, args.refGenome)
                if nlod >= args.nlod_threshold:
                    # hack attempt to satisfy downstream software - not correct for germline heterozygous, male X, etc, etc, etc
                    alts = alt.split(",")
                    if len(alts) == 2:
                        genotype = '1/2'
                    elif len(alts) != 1:
                        raise Exception("error hacking genotype field for " + alts)
                    elif chrom == "chrY" or chrom == "chrM":
                        genotype = '1'
                    elif float(vmf) > 0.95:
                        genotype = '1/1'
                    else:
                        genotype = '0/1'
                    refmt = str(int(umt) - int(vmt))
                    ad = refmt + "," + vmt
                    if len(alts) == 2:
                        ad = ad + ",1"  # horrific hack for the 2nd alt
                    if filter == 'PASS':
                        filter = []
                        if nearby_indel_cnt[int(pos)] >= 3:
                            filter.append('MA')
                        if is_zero_coverage_in_normal:
                            filter.append('ZeroCovInNormal')
                        if is_exist_in_normal:
                            filter.append('ObInNormal')
                        if is_triple_alleles:
                            filter.append('TriAllele')
                        if len(filter) == 0:
                            filter = 'PASS'
                        else:
                            filter = ';'.join(filter)
                    else:
                        if nearby_indel_cnt[int(pos)] >= 3:
                            filter += ';MA'
                        if is_zero_coverage_in_normal:
                            filter += ';ZeroCovInNormal'
                        if is_exist_in_normal:
                            filter += ';ObInNormal'
                        if is_triple_alleles:
                            filter.append('TriAllele')

                    # output
                    format = 'GT:AD:VF'
                    sample = ":".join((genotype, ad, vmf))
                    id = '.'
                    info += ';NLOD=' + str(round(nlod, 2))
                    vcf_line = '\t'.join((chrom, pos, id, ref, alt, qual, filter, info, format, sample)) + '\n'
                    short_line = '\t'.join((chrom, pos, ref, alt, vtype, dp, mt, umt, pi, thr, vmt, vmf, vsm, filter)) + '\n'
                    out_vcf.write(vcf_line)
                    out_variants.write(short_line)

    out_variants.close()
    out_all.close()
    out_variants.close()

    # log run completion
    time_end = datetime.datetime.now()
    print("smCounter completed running at " + str(time_end))
    print("smCounter total time: " + str(time_end - time_start))

    # pass threshold back to caller
    return threshold


# ----------------------------------------------------------------------------------------------
# pythonism to run from the command line
# ----------------------------------------------------------------------------------------------
if __name__ == "__main__":
    # init the argumet parser
    argParseInit()

    # get command line arguments
    args = parser.parse_args()

    # initialize logger
    import run_log

    run_log.init(args.logFile)

    # call main program
    main(args)