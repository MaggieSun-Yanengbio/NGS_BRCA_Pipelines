# run smCounter to variant call
import sys
import os
import time

def sort_index(samfile,samtools):
    bam = samfile[:-3] + 'bam'
    cmd1 = samtools + ' view -b -S ' + samfile + ' > ' + bam
    os.popen(cmd1)
    sorted = bam[:-4] + '_sorted.bam'
    cmd2 = samtools + ' sort ' + bam + ' > ' + sorted
    os.popen(cmd2)
    cmd3 = samtools + ' index ' + sorted
    os.popen(cmd3)
    print ('sort_index has been completed.')
    return sorted

def smCounter(smcounter,outPrefix,bamFile,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
              mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile):
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
    os.popen(cmd)
    print ('smcounter has been completed.')

def main():
    print ('start time is %s.'%time.ctime())
    (source,sample_name,samtools,smcounter,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
     mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath) = sys.argv[1:]
#    add_header(head,samfile)
    outPrefix = source + sample_name + '_variant'
    logFile = source + sample_name + '_logfile'
    samfile = source + 'aligned/' + sample_name + '_vcready.sam'
    bam = sort_index(samfile,samtools)
    smCounter(smcounter,outPrefix,bam,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
              mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile)
    print ('end time is %s.'%time.ctime())

if __name__ == '__main__':
    main()

