# run smCounter to variant call
import sys
import os
import time

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
def call(smcounter,outprefix,bamfile,normal,bedtarget,fastatarget,rpb,ncpu,minbq,minmq,hplen,mismatchthr,
         mtdrop,threshold,min_frequency, min_active_score, tlod_threshold, nlod_threshold,refgenome,
         bedtandemrepeats,bedrepeatmaskersubset,bedtoolspath,logfile):
    cmd = 'python2.7 ' + smcounter + \
          ' --outPrefix ' + outprefix + \
          ' --bamFile ' + bamfile + \
          ' --normal_bamfile ' + normal + \
          ' --bedTarget ' + bedtarget + \
          ' --fastaTarget ' + fastatarget + \
          ' --rpb ' + rpb + \
          ' --nCPU ' + ncpu + \
          ' --minBQ ' + minbq + \
          ' --minMQ ' + minmq + \
          ' --hpLen ' + hplen + \
          ' --mismatchThr ' + mismatchthr + \
          ' --mtDrop ' + mtdrop + \
          ' --threshold ' + threshold + \
          ' --min_frequency ' + min_frequency + \
          ' --min_active_score ' + min_active_score + \
          ' --tlod_threshold ' + tlod_threshold + \
          ' --nlod_threshold ' + nlod_threshold + \
          ' --refGenome ' + refgenome + \
          ' --bedTandemRepeats ' + bedtandemrepeats + \
          ' --bedRepeatMaskerSubset ' + bedrepeatmaskersubset + \
          ' --bedtoolsPath ' + bedtoolspath + \
          ' --logFile ' + logfile
    os.system(cmd)
