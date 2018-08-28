from collections import defaultdict
import pysam
import os
import operator
import random
import scipy.stats
import numpy as np

def get_target_reads(target_region, bamfile, outbam):
    infile = pysam.Samfile(bamfile, 'rb')
    outfile = pysam.Samfile(outbam, 'wb', template=infile)
    for read in infile.fetch(target_region[3:]):
        outfile.write(read)
    outfile.close()


def split_bam(bamfile, half1, half2, split_ratio):
    infile = pysam.Samfile(bamfile, 'rb')
    outfile1 = pysam.Samfile(half1, 'wb', template=infile)
    outfile2 = pysam.Samfile(half2, 'wb', template=infile)
    for read in infile.fetch():
        if random.random() <= split_ratio:
            outfile1.write(read)
        else:
            outfile2.write(read)
    outfile1.close()
    outfile2.close()


def sort_index(bamfile, sorted_bam):
    cmd = 'samtools sort ' + bamfile + ' > ' + sorted_bam
    os.system(cmd)
    cmd = 'samtools index ' + sorted_bam
    os.system(cmd)


def write_fasta(target_region, fasta, refgenome):
    refseq = pysam.FastaFile(refgenome)
    target_seq = refseq.fetch(reference=target_region)
    outfile = open(fasta, 'w')
    header = '>' + target_region + '\n'
    outfile.write(header)
    outfile.write(target_seq)
    outfile.close()


def find_candidate_site(vcf, target_region):
    variants_set = []
    with open(vcf, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                chrom, pos, id, ref, alt = line.split('\t')[:5]
                gt_12878, gt_12889, gt_12890, gt_12891 = line.split('\t')[1340:1344]
                if chrom == target_region and gt_12878 == '0/0' and gt_12891 in ['1/0', '0/1']:
                    chrom = 'chr' + chrom
                    variants_set.append((chrom, pos, ref, alt))
    return (variants_set)


def get_base(read, refseq, chrom, pos):
    if read.indel > 0:
        site = read.alignment.query_sequence[read.query_position]
        inserted = read.alignment.query_sequence[
                   (read.query_position + 1): (read.query_position + 1 + read.indel)]
        base = 'INS|' + site + '|' + site + inserted
    elif read.indel < 0:
        site = read.alignment.query_sequence[read.query_position]
        deleted = refseq.fetch(reference=chrom, start=int(pos), end=int(pos) + abs(read.indel))
        deleted = deleted.upper()
        base = 'DEL|' + site + deleted + '|' + site
    else:
        if read.is_del:
            base = 'DEL'
        else:
            base = read.alignment.query_sequence[read.query_position]
    return base


def filter_variants(variants_set, na12878_bam, na12891_bam, refgenome):
    refseq = pysam.FastaFile(refgenome)
    na12878 = pysam.AlignmentFile(na12878_bam, 'rb')
    na12891 = pysam.AlignmentFile(na12891_bam, 'rb')
    final_variants = []

    for variant in variants_set:
        na12878_base_cnt = defaultdict(int)
        na12891_base_cnt = defaultdict(int)
        na12878_depth_cnt = 0
        na12891_depth_cnt = 0

        chrom, pos, ref, alt = variant
        ref_fa = refseq.fetch(reference=chrom, start=int(pos) - 1, end=int(pos))
        ref_fa = ref_fa.upper()

        na12878_pileup = na12878.pileup(region=chrom[3:] + ':' + pos + ':' + pos, truncate=True, stepper='nofilter')
        na12891_pileup = na12891.pileup(region=chrom[3:] + ':' + pos + ':' + pos, truncate=True, stepper='nofilter')

        for reads_12878 in na12878_pileup:
            for read_12878 in reads_12878.pileups:
                base_12878 = get_base(read_12878, refseq, chrom, pos)
                na12878_base_cnt[base_12878] += 1
                na12878_depth_cnt += 1

        for reads_12891 in na12891_pileup:
            for read_12891 in reads_12891.pileups:
                base_12891 = get_base(read_12891, refseq, chrom, pos)
                na12891_base_cnt[base_12891] += 1
                na12891_depth_cnt += 1

        # filter variants
        if na12878_base_cnt[ref_fa] == na12878_depth_cnt and na12878_depth_cnt >= 200:
            if na12891_depth_cnt > 1 and na12878_depth_cnt >= 100:
                most_common_allele, second_common_allele = sorted(na12891_base_cnt.items(), key=operator.itemgetter(1))[:2]
                if most_common_allele[0] in [ref, alt] and second_common_allele[0] in [ref, alt] and \
                    na12891_base_cnt[second_common_allele[0]]/float(na12891_depth_cnt) >= 0.4 and \
                    (na12891_base_cnt[most_common_allele[0]] + na12891_base_cnt[second_common_allele[0]])/float(na12891_depth_cnt) >= 0.9:
                    final_variants.append(variant)
    return final_variants


def change_reads(final_variants, tumor_half, rate):
    #change reads of filtered variants
    changed_half = pysam.AlignmentFile(tumor_half, 'rb')
    variants_to_depth = defaultdict(int)
    changed_id_to_read = {}

    for variant in final_variants:
        reads_list_half = []
        chrom, pos, ref, alt = variant
        pileup_half = changed_half.pileup(region=chrom[3:] + ':' + pos + ':' + pos, truncate=True, stepper='nofilter')
        for reads_half in pileup_half:
            for read_half in reads_half.pileups:
                reads_list_half.append(read_half)
        variants_to_depth[variant] = len(reads_list_half)

        choosed_reads = random.sample(reads_list_half, int(len(reads_list_half) * rate) + 1)
        for read in reads_list_half:
            if read in choosed_reads:
                qpos = read.query_position
                # if the read has been changed before
                if read.alignment.query_name in changed_id_to_read:
                    read = changed_id_to_read[read.alignment.query_name]
                seq = list(read.alignment.query_sequence)
                qualities = read.alignment.query_qualities
                try:
                    seq[qpos] = alt
                except TypeError:
                    continue
                read.alignment.query_sequence = ''.join(seq)
                read.alignment.query_qualities = qualities
            changed_id_to_read[read.alignment.query_name] = read
    return (changed_id_to_read, variants_to_depth)


def build_final_bam(id_to_read, half1, tumor_half1):
    #write the tumor_half bamfile
    infile = pysam.AlignmentFile(half1, 'rb')
    outfile = pysam.AlignmentFile(tumor_half1, 'wb', template=infile)
    inreads = infile.fetch()
    for read in inreads:
        if read.query_name in id_to_read:
            change_read = id_to_read[read.query_name].alignment
            outfile.write(change_read)
        else:
            outfile.write(read)
    outfile.close()


def recall(smcounter,outprefix,tumor_bam,bed_target,fasta_target,mt_depth,rpb,ncpu,minbq,minmq,hplen,mismatch_thr,mt_drop,
          threshold,refgenome,bed_tandem_repeats,bed_repeat_masker_subset,bedtools_path,logfile):
    #variant calling use the new bamfile
    cmd = 'python2.7 ' + smcounter + \
          ' --outPrefix ' + outprefix + \
          ' --bamFile ' + tumor_bam + \
          ' --bedTarget ' + bed_target + \
          ' --fastaTarget ' + fasta_target + \
          ' --mtDepth ' + mt_depth + \
          ' --rpb ' + rpb + \
          ' --nCPU ' + ncpu + \
          ' --minBQ ' + minbq + \
          ' --minMQ ' + minmq + \
          ' --hpLen ' + hplen + \
          ' --mismatchThr ' + mismatch_thr + \
          ' --mtDrop ' + mt_drop + \
          ' --threshold ' + threshold + \
          ' --refGenome ' + refgenome + \
          ' --bedTandemRepeats ' + bed_tandem_repeats + \
          ' --bedRepeatMaskerSubset ' + bed_repeat_masker_subset + \
          ' --bedtoolsPath ' + bedtools_path + \
          ' --logFile ' + logfile
    os.system(cmd)


def check_half2(half2, chrom, pos, ref, alt, refgenome):
    # check if has alt base in normal sample bamfile at FP sites.
    refseq = pysam.FastaFile(refgenome)
    alt_cnt = 0
    depth_cnt = 0
    normal = pysam.AlignmentFile(half2, 'rb')
    if len(ref) > 1 and len(alt) == 1:
        for reads in normal.pileup(region=chrom[3:] + ':' + pos + ':' + pos, truncate=True, stepper='nofilter'):
            for read in reads.pileups:
                depth_cnt += 1
                if read.indel < 0:
                    deleted = refseq.fetch(reference=chrom, start=int(pos), end=int(pos) + abs(read.indel))
                    if deleted == ref[1:]:
                        alt_cnt += 1
    elif len(ref) == 1 and len(alt) > 1:
        for reads in normal.pileup(region=chrom[3:] + ':' + pos + ':' + pos, truncate=True, stepper='nofilter'):
            for read in reads.pileups:
                depth_cnt += 1
                if read.indel > 0:
                    inserted = read.alignment.query_sequence[(read.query_position + 1): (read.query_position + 1 + read.indel)]
                    if inserted == alt[1:]:
                        alt_cnt += 1
    else:
        for reads in normal.pileup(region=chrom[3:] + ':' + pos + ':' + pos, truncate=True, stepper='nofilter'):
            for read in reads.pileups:
                depth_cnt += 1
                if read.indel == 0:
                    base = read.alignment.query_sequence[read.query_position]
                    if base == alt:
                        alt_cnt += 1

    if depth_cnt == 0:
        return 'Zero_Depth'
    else:
        return alt_cnt


def read_new_variants(half2, refgenome, vcf, variants_list, detail_table, rate, sample, depth_true, cycle, is_count_fp, is_check_tp):
    #statistics new vcf
    tp_cnt = 0
    pass_cnt = 0
    pass_variants = []
    for l1 in open(vcf, 'r'):
        if not l1.startswith('#'):
            chrom, pos, xxx, ref, alt, qual, flt, info, fmt, freq = l1.split('\t')
            if flt == 'PASS':
                tlod, nlod = info.strip().split(';')[-2:]
                tlod = tlod[tlod.find('=')+1:]
                nlod = nlod[nlod.find('=')+1:]
                gt, dp, fr = freq.strip().split(':')
                dp = dp.replace(',',';')
                pass_cnt += 1
                alt_cnt = check_half2(half2, chrom, pos, ref, alt, refgenome)
                pass_variants.append((str(cycle+1), chrom, pos, ref, alt, dp, tlod, nlod, fr, str(alt_cnt)))

                if (chrom, pos, ref, alt) in variants_list:
                    if is_check_tp:
                        detail_table.write(','.join([sample, str(rate), str(cycle+1), chrom, pos, ref, alt, dp, tlod, nlod ,fr]) + '\n')
                    ref_dp, alt_dp = dp.split(';')
                    depth_true[(chrom, pos, ref, alt)]['ref'].append(int(ref_dp))
                    depth_true[(chrom, pos, ref, alt)]['alt'].append(int(alt_dp))
                    depth_true[(chrom, pos, ref, alt)]['freq'].append(float(fr))
                    depth_true[(chrom, pos, ref, alt)]['tlod'].append(float(tlod))
                    depth_true[(chrom, pos, ref, alt)]['nlod'].append(float(nlod))
                    tp_cnt += 1
    if is_count_fp:
        return (pass_cnt, pass_variants)
    if is_check_tp:
        return (depth_true, tp_cnt)


def main():
    path = '/home/administrator/benchmark_test/'
    sample = 'NA12878'
    target_region = 'chr17'
    caller = '/home/administrator/smCounter/smCounter-master/caller_na12878.py'
    bed_target = "/home/administrator/source/target_breast.bed"
    mt_depth = '300'
    rpb = '8'
    ncpu = '16'
    minbq = '20'
    minmq = '20'
    hplen = '8'
    mismatch_thr = '6'
    mtdrop = '0'
    threshold = '7'
    refgenome = "/home/administrator/genome/ucsc.hg19.fasta"
    bed_tandem_repeats = "/home/administrator/smCounter/smCounter-master/simpleRepeat.bed"
    bed_repeat_masker_subset = "/home/administrator/smCounter/smCounter-master/SR_LC_SL.nochr.bed"
    bedtools_path = "/usr/bin/"
    logfile = '/home/administrator/benchmark_test/na12878_logfile'

    tp_summary = open(path + sample + '_' + target_region + '_tp_summary.csv', 'w')
    tp_detail_table = open(path + sample + '_' + target_region + '_tp_detail_table.csv', 'w')
    tp_called_count = open(path + sample + '_' + target_region + '_tp_called_count.csv', 'w')
    fp_total_variants = open(path + sample + '_' + target_region + '_fp_total_variants.csv', 'w')
    fp_per_cycle = open(path + sample + '_' + target_region + '_fp_per_cycle.csv', 'w')
    fp_variants_count = open(path + sample + '_' + target_region + '_fp_variants_count.csv', 'w')

    tp_summary.write(','.join(['sample','frequency_theory','chromosome','position','ref','alt','tlod','nlod','true_ref_depth','true_alt_depth','frequency_true','pvalue\n']))
    tp_detail_table.write(','.join(['sample','frequency_theory','cycle','chromosome','position','ref','alt','depth','tlod','nlod','frequency_true\n']))
    tp_called_count.write(','.join(['frequency_theory', 'tp_true', 'tp_theory', 'ratio\n']))
    fp_total_variants.write(','.join(['cycle', 'chromosome', 'position', 'ref', 'alt', 'depth', 'tlod', 'nlod', 'frequency', 'exist_alt_in_normal\n']))
    fp_per_cycle.write('exp_idx,fp_cnt\n')
    fp_variants_count.write(','.join(['chromosome', 'position', 'ref', 'alt', 'count\n']))

    na12878_bam = '/home/administrator/benchmark_test/NA12878_' + target_region + '.bam'
    na12878_sorted = na12878_bam[:-4] + '_sorted.bam'
    na12891_bam = '/home/administrator/benchmark_test/NA12891_' + target_region + '.bam'
    na12891_sorted = na12891_bam[:-4] + '_sorted.bam'
    na12878_half1 = na12878_bam[:-4] + '_half1.bam'
    na12878_half1_sorted = na12878_half1[:-4] + '_sorted.bam'
    na12878_half2 = na12878_bam[:-4] + '_half2.bam'
    na12878_half2_sorted = na12878_half2[:-4] + '_sorted.bam'
    na12878_tumor_half1 = na12878_bam[:-4] + '_tumor_half1.bam'

    #extrace target data in fasta, bam and vcf
    fasta_target = path + target_region + '.fa'
    write_fasta(target_region, fasta_target, refgenome)
    get_target_reads(target_region, '/home/administrator/benchmark_test/NA12878.mapped.illumina.mosaik.CEU.exome.20110411.bam', na12878_bam)
    sort_index(na12878_bam, na12878_sorted)
    get_target_reads(target_region, '/home/administrator/benchmark_test/NA12891.mapped.illumina.mosaik.CEU.exome.20110411.bam', na12891_bam)
    sort_index(na12891_bam, na12891_sorted)

    variants_set = find_candidate_site('/home/administrator/benchmark_test/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf', target_region[3:])
    fp_total_cnt = defaultdict(int)
    fp_cnt_list = []
    exp_cnt = 0

    for rate in [0.01, 0.02, 0.03, 0.04, 0.05]:
        outprefix_fp = '/home/administrator/benchmark_test/NA12878_half1_fp'
        outprefix_tp = '/home/administrator/benchmark_test/NA12878_half1_tp_' + str(rate)
        na12878_tumor_half1_sorted = na12878_tumor_half1[:-4] + '_' + str(rate) + '_sorted.bam'
        #filter variants by depth, frequency and base
        final_variants = filter_variants(variants_set, na12878_sorted, na12891_sorted,  refgenome)
        tp_cnt_list = []
        depth_true = defaultdict(lambda: defaultdict(list))
        cycle = 0

        while cycle < 10:
            # divide na12878 data into 2 parts.
            split_bam(na12878_sorted, na12878_half1, na12878_half2, 0.5)
            sort_index(na12878_half1, na12878_half1_sorted)
            sort_index(na12878_half2, na12878_half2_sorted)

            # variant calling, get false positive.
            recall(caller, outprefix_fp, na12878_half1_sorted, bed_target, fasta_target, mt_depth, rpb, ncpu, minbq,
                   minmq, hplen, mismatch_thr, mtdrop, threshold, refgenome, bed_tandem_repeats, bed_repeat_masker_subset,
                   bedtools_path, logfile)

            # variant calling, get virtual somatic mutations.
            changed_id_to_reads, variants_to_depth = change_reads(final_variants, na12878_half1_sorted, rate)
            build_final_bam(changed_id_to_reads, na12878_half1_sorted, na12878_tumor_half1)
            sort_index(na12878_tumor_half1, na12878_tumor_half1_sorted)
            recall(caller, outprefix_tp, na12878_tumor_half1_sorted, bed_target, fasta_target, mt_depth, rpb, ncpu,
                   minbq, minmq, hplen, mismatch_thr, mtdrop, threshold, refgenome, bed_tandem_repeats, bed_repeat_masker_subset,
                   bedtools_path, logfile)

            depth_true, tp_cnt = read_new_variants(na12878_half2_sorted, refgenome, outprefix_tp+'.smCounter.somatic.cut.vcf', final_variants,
                                                   tp_detail_table, rate, sample, depth_true, cycle, False, True)
            fp_cnt, fp_variants = read_new_variants(na12878_half2_sorted, refgenome, outprefix_fp + '.smCounter.somatic.cut.vcf', final_variants,
                                                    tp_detail_table, rate, sample, depth_true, cycle, True, False)
            for fp_variant in fp_variants:
                fp_total_variants.write(','.join(fp_variant) + '\n')
                fp_total_cnt[fp_variant[1:5]] += 1
            fp_cnt_list.append(fp_cnt)
            tp_cnt_list.append(tp_cnt)

            exp_cnt += 1
            fp_per_cycle.write(','.join([str(exp_cnt), str(fp_cnt)]) + '\n')
            cycle += 1
            print('***********************************************************')
            print('Rate:%s, Cycle %s finish.' % (rate, cycle))
            print('***********************************************************')

        tp_mean = np.mean(tp_cnt_list)
        tp_called_count.write(','.join([str(rate), str(tp_mean), str(len(final_variants)), str(tp_mean/len(final_variants))]) + '\n')

        for variant in sorted(depth_true):
            options_to_list = depth_true[variant]
            chrom, pos, ref, alt = variant
            frequency = np.mean(options_to_list['freq'])
            true_ref = np.mean(options_to_list['ref'])
            true_alt = np.mean(options_to_list['alt'])
            tlod = np.mean(options_to_list['tlod'])
            nlod = np.mean(options_to_list['nlod'])
            theory_ref = int(variants_to_depth[variant] * (1 - rate))
            theory_alt = int(variants_to_depth[variant] * rate)
            oddsratio, pvalue = scipy.stats.fisher_exact([[theory_ref, theory_alt], [int(true_ref), int(true_alt)]])
            tp_summary.write(','.join([sample, str(rate), chrom, pos, ref, alt, str(tlod), str(nlod), str(true_ref),
                                str(true_alt), str(frequency), str(pvalue) +'\n']))

    for fp_vrt, cnt in sorted(fp_total_cnt.items(),key = operator.itemgetter(1)):
        chrom, pos, ref, alt = fp_vrt
        fp_variants_count.write(','.join([chrom, pos, ref, alt, str(cnt)]) + '\n')
    fp_per_cycle.write(','.join(['mean', str(np.mean(fp_cnt_list))]) + '\n')

    tp_detail_table.close()
    tp_summary.close()
    tp_called_count.close()
    fp_per_cycle.close()
    fp_variants_count.close()
    fp_total_variants.close()

if __name__ == '__main__':
    main()
