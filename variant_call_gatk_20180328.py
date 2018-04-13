##somatic variant call by gatk-mutect2
##input: samfile(aligned)
##output: variant.vcf
import os
import sys

def sort_markdups(sample,samfile,sorted,bam,markdups,metrics,add_RG,gatk,samtools,genome):
    cmd1 = 'java -jar '+gatk+' SortSam -I '+samfile+' -O '+sorted+' -SO coordinate'
    os.system(cmd1)
    # Prepping reference geonome fasta
    cmd2 = 'java -jar '+gatk+' CreateSequenceDictionary -R '+genome
    os.system(cmd2)
    cmd3 = samtools+' faidx '+genome
    os.system(cmd3)
    cmd4 = samtools+' view -bS '+sorted+' > '+bam
    os.system(cmd4)
    cmd5 = 'java -jar '+gatk+' MarkDuplicates -I '+bam+' -O '+markdups+' -M '+metrics+' --REMOVE_SEQUENCING_DUPLICATES false'
    os.system(cmd5)
    #add row 'RG' in head
    cmd6 = 'java -jar '+gatk+' AddOrReplaceReadGroups -I '+markdups+' -O '+add_RG+' -LB lib1 -PL illumina -PU unit1 -SM '+sample
    os.system(cmd6)
    cmd7 = samtools+' index '+add_RG
    os.system(cmd7)
    return add_RG

def IndexFeatureFile(knownsites,gatk):
    cmd1 = 'java -jar '+gatk+' IndexFeatureFile -F '+knownsites
    os.system(cmd1)

def ApplyBQSR(bamfile,recal_table,recal_bam,gatk,genome,knownsites):
    cmd1 = 'java -jar '+gatk+' BaseRecalibrator -I '+bamfile+' -O '+recal_table+' -R '+genome+knownsites
    os.system(cmd1)
    cmd2 = 'java -jar '+gatk+' ApplyBQSR -I '+bamfile+' -O '+recal_bam+' -bqsr '+recal_table
    os.system(cmd2)
    return recal_bam

def Mutect2_call(recal_bam,raw_variant,variant_filtered,gatk,genome,sample):
    cmd1 = 'java -jar '+gatk+' Mutect2 -I '+recal_bam+' -O '+raw_variant+' -R '+genome+' -tumor '+sample
    os.system(cmd1)
    cmd2 = 'java -jar '+gatk+' FilterMutectCalls -V '+raw_variant+' -O '+variant_filtered
    os.system(cmd2)

def main():
    source, sample, gatk, samtools, genome = sys.argv[1:6]
    knownsites = sys.argv[6:]
    source += 'aligned/'
    samfile = source + sample + '_formated.sam'
    sorted = source + sample + '_sorted.sam'
    bam = source + sample + '_sorted.bam'
    markdups = source + sample + '_dups_marked.bam'
    metrics = source + sample + '_dups_metrics.txt'
    add_RG = source + sample + '_dups_marked_RG.bam'
    recal_ready = sort_markdups(sample,samfile,sorted,bam,markdups,metrics,add_RG,gatk,samtools,genome)
    cmd_knownsites = ''
    for site in knownsites:
        IndexFeatureFile(site, gatk)
        cmd_knownsites += ' --known-sites '+site+' '
    recal_table = source + sample + '_recal_data.table'
    recal_bam = source + sample + '_recal_reads.bam'
    vc_ready = ApplyBQSR(recal_ready, recal_table, recal_bam, gatk, genome, cmd_knownsites)
    raw_variant = source + sample + '_raw_variant.vcf'
    variant_filtered = source + sample + '_variant_filtered.vcf'
    Mutect2_call(vc_ready, raw_variant, variant_filtered, gatk, genome, sample)

if __name__ == '__main__':
    main()
