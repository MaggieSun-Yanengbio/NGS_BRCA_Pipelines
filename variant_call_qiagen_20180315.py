# run smCounter to variant call
import sys
import os
import time

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
    print ('sort_index has been completed.')
    return sorted

def smCounter(smcounter,outPrefix,bamFile,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
              mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile):
    #outPrefix--prefix for output files
    #bamFile--Bam file
    #bedTarget--BED file for target region
    #mtDepth--Mean MT depth
    #rpb--Mean read pairs per MT
    #nCPU--number of CPUs to use in parallel
    #minBQ--minimum base quality allowed for analysis
    #minMQ--minimum mapping quality allowed for analysis
    #hpLen--Minimum length of homopolymers
    #mismatchThr--average number of mismatches per 100 bases allowed
    #mtDrop--Drop MTs with lower than or equal to X reads.
    #threshold--Minimum prediction index for a variant to be called.Must be non-negative. Typically ranges from 10 to 60.If set to 0 (default), smCounter will choose the appropriate cutoff based on the mean MT depth.
    #refGenome--refgenome
    #bedTandemRepeats--bed for UCSC tandem repeats
    #bedRepeatMaskerSubset--bed for RepeatMasker simple repeats, low complexity,microsatellite regions
    #bedtoolsPath--path to bedtools
    #logFile--log file
    cmd = 'python2.7 ' + smcounter + \
          ' --outPrefix ' + outPrefix + \
          ' --bamFile ' + bamFile + \
          ' --bedTarget ' +bedTarget + \
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
    os.system(cmd)
    print ('smcounter has been completed.')

def main():
    time_start = time.time()
    (source,sample_name,samtools,smcounter,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
     mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath) = sys.argv[1:]
#    add_header(head,samfile)
    outPrefix = source + sample_name + '_variant'
    logFile = source + sample_name + '_logfile'
    samfile = source + 'aligned/' + sample_name + '_vcready.sam'
    bam = sam_to_bam(samfile,samtools)
    sorted = sort_index(bam,samtools)
    smCounter(smcounter,outPrefix,sorted,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
              mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile)
    print ('The time of used smcounter is %s minutes.'%str((time.time()-time_start)/60))

if __name__ == '__main__':
    main()

