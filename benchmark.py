import os
import pysam
from collections import defaultdict
import random
import numpy as np
import scipy.stats

def readVcf(smcounter_vcf, mutect2_vcf, varscan2_vcf, frequency):
    # read the raw vcf,get variants: 1.frequency>=95%
    #                                2.depth>=10
    #                                3.these variants in smcounter.vcf mutect2.vcf and varscan2.vcf
    depth_theory = defaultdict(lambda: defaultdict(int))
    list1 = []
    list2 = []
    list3 = []
    for l1 in open(smcounter_vcf, 'r'):
        if not l1.startswith('#'):
            chrom, pos, xxx, ref, alt, qual, fil, info, fmt, freq = l1.split('\t')
            gt, dp, fr = freq.split(':')
            ref_dp, alt_dp = dp.split(',')
            if float(fr) >= frequency:
                list1.append((chrom, pos, ref, alt))
                depth_theory[(chrom, pos, ref, alt)]['smcounter'] = int(ref_dp) + int(alt_dp)

    for l2 in open(mutect2_vcf, 'r'):
        if not l2.startswith('#'):
            chrom, pos, xxx, ref, alt, qual, fil, info, fmt, freq = l2.split('\t')
            dp = int(info[info.find('DP')+3:info.find(';')])
            if ',' in alt:
                alt = alt.split(',')[0]
                if float(freq.split(':')[2].split(',')[0]) >= frequency and dp >= 10:
                    list2.append((chrom, pos, ref ,alt))
                    depth_theory[(chrom, pos, ref, alt)]['mutect2'] = int(dp)
            else:
                if float(freq.split(':')[2]) >= frequency and dp >= 10:
                    list2.append((chrom, pos, ref, alt))
                    depth_theory[(chrom, pos, ref, alt)]['mutect2'] = int(dp)

    f3 = open(varscan2_vcf, 'r')
    f3.readline()
    for l3 in f3.readlines():
        chrom, pos, ref, cons, r1, r2, freq = l3.split('\t')[:7]
        alt = l3.strip().split('\t')[-1]
        if float(freq[:-1])/100 >= frequency:
            list3.append((chrom, pos, ref, alt))
            depth_theory[(chrom, pos, ref, alt)]['varscan2'] = int(r1) + int(r2)

    variants_list = [variant for variant in list1 if variant in list2 and variant in list3]

    return variants_list, depth_theory

def changeToReference(infile, variants_list, rate):
    #randomly choose reads*rate
    #change to reference base
    id_to_read = {}
    for variant in variants_list:
        reads = []
        chrom, pos, ref, alt = variant
        inreads = infile.fetch(chrom, int(pos)-1, int(pos))
        for read in inreads:
            reads.append(read)
        choosed_reads = random.sample(reads, int(len(reads)*rate))
        for read in choosed_reads:
            seq = list(read.query_sequence)
            qualities = read.query_qualities
            if read.cigartuples[0][0] == 4:
                seq[int(pos)-read.pos+read.cigartuples[0][1]-1] = ref
            else:
                seq[int(pos)-read.pos-1] = ref
            read.query_sequence = ''.join(seq)
            read.query_qualities = qualities
            id_to_read[read.query_name] = read

    return id_to_read

def buildFinalBam(id_to_read, infile, outfile):
    #write the changed bamfile
    inreads = infile.fetch()
    for read in inreads:
        if read.query_name in id_to_read:
            read = id_to_read[read.query_name]
        outfile.write(read)
    outfile.close()

def reCall(smcounter,outPrefix,changed_bam,bedTarget,fastaTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,mtDrop,
          threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile,
          gatk,gatk_sorted,add_RG,knownsites1,knownsites2,knownsites3,recal_table,recal_bam,raw_variant,variant_filtered,
          samtools, varscan, source, sample, change):
    #variant calling use the new bamfile
    #smcounter
    cmd1 = 'python2.7 ' + smcounter + \
          ' --outPrefix ' + outPrefix + \
          ' --bamFile ' + changed_bam + \
          ' --bedTarget ' + bedTarget + \
          ' --fastaTarget ' + fastaTarget + \
          ' --mtDepth ' + mtDepth + \
          ' --rpb ' + rpb + \
          ' --nCPU ' + nCPU + \
          ' --minBQ ' + minBQ + \
          ' --minMQ ' + minMQ + \
          ' --hpLen ' + hpLen + \
          ' --mismatchThr ' + mismatchThr + \
          ' --mtDrop ' + mtDrop + \
          ' --threshold ' + threshold + \
          ' --refGenome ' + refGenome + \
          ' --bedTandemRepeats ' + bedTandemRepeats + \
          ' --bedRepeatMaskerSubset ' + bedRepeatMaskerSubset + \
          ' --bedtoolsPath ' + bedtoolsPath + \
          ' --logFile ' + logFile
    os.system(cmd1)

    #mutect2
    cmd2 = 'java -jar ' + gatk + ' SortSam -I ' + changed_bam + ' -O ' + gatk_sorted + ' -SO coordinate'
    os.system(cmd2)
    cmd3 = 'java -jar ' + gatk + ' AddOrReplaceReadGroups -I ' + gatk_sorted + ' -O ' + add_RG + ' -LB lib1 -PL illumina -PU unit1 -SM ' + sample
    os.system(cmd3)
    cmd4 = samtools + ' index ' + add_RG
    os.system(cmd4)
    cmd5 = 'java -jar ' + gatk + ' BaseRecalibrator -I ' + add_RG + ' -O ' + recal_table + ' -R ' + refGenome + ' --known-sites ' + knownsites1 + ' --known-sites ' + knownsites2 + ' --known-sites ' + knownsites3
    os.system(cmd5)
    cmd6 = 'java -jar ' + gatk + ' ApplyBQSR -I ' + add_RG + ' -O ' + recal_bam + ' -bqsr ' + recal_table
    os.system(cmd6)
    cmd7 = 'java -jar ' + gatk + ' Mutect2 -I ' + recal_bam + ' -O ' + raw_variant + ' -R ' + refGenome + ' -tumor ' + sample
    os.system(cmd7)
    cmd8 = 'java -jar ' + gatk + ' FilterMutectCalls -V ' + raw_variant + ' -O ' + variant_filtered
    os.system(cmd8)

    #varscan2
    cmd9 = samtools + ' mpileup -f ' + refGenome + ' ' + changed_bam + ' > ' + source + sample + '_pileup'
    os.system(cmd9)
    cmd10 = 'java -jar ' + varscan + ' pileup2snp ' + source + sample + '_pileup ' + '--min-coverage 8 --min-reads2 5 > ' + source + sample + '_snp_' + change
    os.system(cmd10)

def getNewVariants(smcounter_vcf, mutect2_vcf, varscan2_vcf, variants_list, detail_table, rate, sample, depth_true, cycle):
    #statistics new vcf
    for l1 in open(smcounter_vcf, 'r'):
        if not l1.startswith('#'):
            chrom, pos, xxx, ref, alt, qual, fil, info, fmt, freq = l1.split('\t')
            gt, dp, fr = freq.split(':')
            dp = dp.replace(',',';')
            if (chrom, pos, ref, alt) in variants_list:
                detail_table.write(','.join([sample, 'smcounter', str(rate), str(cycle+1), chrom, pos, ref, alt, dp, fr]))
                ref_dp, alt_dp = dp.split(';')
                depth_true[(chrom,pos,ref,alt)]['smcounter']['ref'].append(int(ref_dp))
                depth_true[(chrom,pos,ref,alt)]['smcounter']['alt'].append(int(alt_dp))
                depth_true[(chrom,pos,ref,alt)]['smcounter']['freq'].append(float(fr))

    for l2 in open(mutect2_vcf, 'r'):
        if not l2.startswith('#'):
            chrom, pos, xxx, ref, alt, qual, fil, info, fmt, freq = l2.split('\t')
            dp = freq.split(':')[1].replace(',',';')
            fr = freq.split(':')[2]
            if (chrom, pos, ref, alt) in variants_list:
                detail_table.write(','.join([sample, 'mutect2', str(rate), str(cycle+1), chrom, pos, ref, alt, dp, fr+'\n']))
                ref_dp, alt_dp = dp.split(';')
                depth_true[(chrom,pos,ref,alt)]['mutect2']['ref'].append(int(ref_dp))
                depth_true[(chrom,pos,ref,alt)]['mutect2']['alt'].append(int(alt_dp))
                depth_true[(chrom,pos,ref,alt)]['mutect2']['freq'].append(float(fr))

    f3 = open(varscan2_vcf, 'r')
    f3.readline()
    for l3 in f3.readlines():
        chrom, pos, ref, cons, r1, r2, freq = l3.split('\t')[:7]
        alt = l3.strip().split('\t')[-1]
        dp = r1 + ';' + r2
        fr = str(float(freq[:-1])/100)
        if (chrom, pos, ref, alt) in variants_list:
            detail_table.write(','.join([sample, 'varscan2', str(rate), str(cycle+1), chrom, pos, ref, alt, dp, fr+'\n']))
            ref_dp, alt_dp = dp.split(';')
            depth_true[(chrom,pos,ref,alt)]['varscan2']['ref'].append(int(ref_dp))
            depth_true[(chrom,pos,ref,alt)]['varscan2']['alt'].append(int(alt_dp))
            depth_true[(chrom,pos,ref,alt)]['varscan2']['freq'].append(float(fr))

    return depth_true

def main():
    path = '/home/administrator/benchmark_test/'
    samples = ['N223-M2_S10', 'N223-M3_S11', 'N223-N_S9', 'Q223-M2_S7', 'Q223-M3_S8', 'Q5000-M2_S5', 'Q5000-M3_S6', 'Q5000-N_S4', 'QIA-M2_S2', 'QIA-M3_S3', 'QIA-N_S1']
    detail_table = open(path + 'detail_table.csv', 'w')
    stats_table = open(path + 'stats.csv', 'w')
    for sample in samples:
        bam = path + sample + '_vcready_sorted.bam'
        smcounter_vcf = path + sample + '_variant.smCounter.cut.vcf'
        mutect2_vcf = path + sample + '_variant_filtered.vcf'
        varscan2_snp = path + sample + '_snp'
        #smcounter parameters
        smcounter = '/home/administrator/smCounter/smCounter-master/smCounter6.py'
        bedTarget = "/home/administrator/source/target_breast.bed"
        fastaTarget = '/home/administrator/source/target_breast.refSeq.fa'
        mtDepth = '300'
        rpb = '8'
        nCPU = '16'
        minBQ = '20'
        minMQ = '30'
        hpLen = '8'
        mismatchThr = '6'
        mtDrop = '0'
        threshold = '7'
        refGenome = "/home/administrator/genome/ucsc.hg19.fasta"
        bedTandemRepeats = "/home/administrator/smCounter/smCounter-master/simpleRepeat.bed"
        bedRepeatMaskerSubset = "/home/administrator/smCounter/smCounter-master/SR_LC_SL.nochr.bed"
        bedtoolsPath = "/usr/bin/"
        logFile = path + sample + '_logfile'
        #mutect2 parameters
        gatk = '/home/administrator/gatk/gatk-package-4.0.1.2-local.jar'
        knownsites1 = '/home/administrator/knowkites/1000G_phase1.indels.hg19.vcf'
        knownsites2 = '/home/administrator/knowkites/1000G_phase1.snps.high_confidence.hg19.vcf'
        knownsites3 = '/home/administrator/knowkites/Mills_and_1000G_gold_standard.indels.hg19.vcf'
        #varscan2 parameters
        samtools = 'samtools'
        varscan = '/home/administrator/varscan/VarScan.v2.4.3.jar'

        infile = pysam.Samfile(bam, 'rb')
        variants_list, depth_theory = readVcf(smcounter_vcf, mutect2_vcf, varscan2_snp, 0.95)

        for rate in [0.99, 0.95, 0.85, 0.65, 0.5]:

            depth_true = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            changed_bam = path + sample + '_change' + str(rate) + '.bam'
            outPrefix = path + sample + '_' + str(rate)
            gatk_sorted = path + sample + '_' + str(rate) + '_gatk_sorted.bam'
            add_RG = path + sample + '_' + str(rate) + '_dups_marked_RG.bam'
            recal_table = path + sample + '_' + str(rate) + '_recal_data.table'
            recal_bam = path + sample + '_' + str(rate) + '_recal_reads.bam'
            raw_variant = path + sample + '_' + str(rate) + '_raw_variant.vcf'
            variant_filtered = path + sample + '_' + str(rate) + '_variant_filtered.vcf'

            cycle = 0
            while cycle < 10:
                outfile = pysam.Samfile(changed_bam, 'wb', template=infile)
                id_to_read = changeToReference(infile, variants_list, rate)
                buildFinalBam(id_to_read, infile, outfile)
                sorted_bam = changed_bam[:-4] + '_sorted.bam'
                cmd1 = samtools + ' sort ' + changed_bam + ' > ' + sorted_bam
                os.system(cmd1)
                cmd2 = samtools + ' index ' + sorted_bam
                os.system(cmd2)
                reCall(smcounter,outPrefix,sorted_bam,bedTarget,fastaTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile,
                       gatk, gatk_sorted, add_RG, knownsites1, knownsites2, knownsites3, recal_table, recal_bam, raw_variant,variant_filtered,
                       samtools, varscan, path, sample, str(rate))
                depth_true = getNewVariants(outPrefix+'.smCounter.cut.vcf', variant_filtered, path+sample+'_snp_'+str(rate), variants_list, detail_table, rate, sample, depth_true, cycle)
                cycle += 1
                print ('***********************************************************')
                print ('Sample:%s, Rate:%s\nCycle %s finish.'%(sample, rate, cycle))
                print ('***********************************************************')

            for variant, software_to_option in depth_true.items():
                chrom, pos, ref, alt = variant
                for software, option in software_to_option.items():
                    frequency = np.mean(option['freq'])
                    true_ref = np.mean(option['ref'])
                    true_alt = np.mean(option['alt'])
                    theory_ref = int(depth_theory[variant][software]*(1-rate))
                    theory_alt = int(depth_theory[variant][software]*rate)
                    oddsratio, pvalue = scipy.stats.fisher_exact([[theory_ref, int(true_ref)], [theory_alt, int(true_alt)]])
                    stats_table.write(','.join([sample, software, str(rate), chrom, pos, ref, alt, str(true_ref), str(true_alt), str(frequency), str(pvalue)+'\n']))
    detail_table.close()
    stats_table.close()

if __name__ == '__main__':
    main()