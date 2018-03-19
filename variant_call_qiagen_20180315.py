# run smCounter to variant call
import sys
import os

def add_header(head, samfile):
    cmd1 = 'cat ' + samfile + ' >> ' + head
    os.popen(cmd1)
    cmd2 = 'rm ' + samfile
    os.popen(cmd2)
    cmd3 = 'mv ' + head + ' ' + samfile
    os.popen(cmd3)
    print ('add_header is complete.')

def sort_index(samfile,samtools):
    bam = samfile[:-3] + 'bam'
    cmd1 = samtools + ' view -b -S ' + samfile + ' > ' + bam
    os.popen(cmd1)
    sorted = bam[:-4] + '_sorted.bam'
    cmd2 = samtools + ' sort ' + bam + ' > ' + sorted
    os.popen(cmd2)
    cmd3 = samtools + ' index ' + sorted
    os.popen(cmd3)
    print ('sort_index is complete.')
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
    print ('smcounter is complete.')

def main():
    (head,samfile,samtools,smcounter,outPrefix,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
     mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile) = sys.argv[1:]
    add_header(head,samfile)
    bam = sort_index(samfile,samtools)
    smCounter(smcounter,outPrefix,bam,bedTarget,mtDepth,rpb,nCPU,minBQ,minMQ,hpLen,mismatchThr,
              mtDrop,threshold,refGenome,bedTandemRepeats,bedRepeatMaskerSubset,bedtoolsPath,logFile)

if __name__ == '__main__':
    main()
