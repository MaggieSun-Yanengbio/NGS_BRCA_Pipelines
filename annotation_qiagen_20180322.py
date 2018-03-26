import sys
import time

#INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">
#INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
#INFO=<ID=MT,Number=1,Type=Integer,Description="Total MT depth">
#INFO=<ID=UMT,Number=1,Type=Integer,Description="Filtered MT depth">
#INFO=<ID=PI,Number=1,Type=Float,Description="Variant prediction index">
#INFO=<ID=THR,Number=1,Type=Integer,Description="Variant prediction index minimum threshold">
#INFO=<ID=VMT,Number=1,Type=Integer,Description="Variant MT depth">
#INFO=<ID=VMF,Number=1,Type=Float,Description="Variant MT fraction">
#INFO=<ID=VSM,Number=1,Type=Integer,Description="Variant strong MT depth">
#FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">
#FORMAT=<ID=VF,Number=1,Type=Float,Description="Variant MT fraction, same as VMF">
#CLNDN--ClinVar`s preferred disease name for the concept specified by disease identifiers in CLNDISD.
#HGVS--Top-level(primary assembly,alt,or patch) HGVS expression.
#CLNSIG--Clinical significance for this single variant.
#Mutation_Description--Type of mutation at the amino acid level (substitution, deletion, insertion, complex, fusion, unknown etc.)
#Gene_CDS_Length--Length of the gene (base pair) counts.
#Mutation_Zygosity--Information on whether the mutation was reported to be homozygous , heterozygous or unknown within the sample.
#LOH--LOH Information on whether the gene was reported to have loss of heterozygosity in the sample: yes, no or unknown.
#Mutation_Strand--postive or negative.
#FATHMM_Prediction--Functional Analysis through Hidden Markov Models.
#FATHMM_Score--The scores are in the form of pvalues ranging from 0 to 1. Pvalues greater than 0.5 are pathogenic
             # while less than 0.5 are benign. Pvalues close to 0 or 1 are the high confidence results which
             # are more accurate. The results are annotated as 10 feature groups (separately for coding and
             #  non coding variants) details of which can be found in the original FATHMM-MKL paper.
#Mutation_Somatic_Status--Information on whether the sample was reported to be Confirmed Somatic, Previously Reported or Variant of unknown origin.

def read_database(cosmic,clinvar):
    cos = open(cosmic, 'r')
    dict_cos = {}
    cos.readline()
    for db1 in cos.readlines():
        # if not db1.startswith('Gene_name'):
        Gene_name, Accession_Number, Gene_CDS_length, HGNC_ID, Sample_name, ID_sample, ID_tumour, Primary_site, \
        Site_subtype1, Site_subtype2, Site_subtype3, Primary_histology, Histology_subtype1, Histology_subtype2, \
        Histology_subtype3, Genome_wide_screen, Mutation_ID, Mutation_CDS, Mutation_AA, Mutation_Description, \
        Mutation_zygosity, LOH, GRCh, Chr, Start, End, Mutation_strand, SNP, Resistance_Mutation, FATHMM_prediction, \
        FATHMM_score, Mutation_somatic_status, Pubmed_PMID, ID_STUDY, Sample_source, Tumor_origin, Age = db1.strip().split(',')
        key1 = [Chr, Start]
        value1 = [HGNC_ID, Mutation_ID, Mutation_Description, Accession_Number, Gene_name, Gene_CDS_length,
                  Mutation_zygosity, LOH, Mutation_strand, Mutation_CDS, Mutation_AA, FATHMM_prediction, FATHMM_score,
                  Mutation_somatic_status]
        dict_cos[','.join(key1)] = ','.join(value1)
    clin = open(clinvar, 'r')
    dict_clin = {}
    clin.readline()
    for db2 in clin.readlines():
        genename, chr, start, end, geneid, alleleid, rs, pos, reference_seq, variant_seq, af_esp, af_exac, af_tgp, \
        clndn, clnhgvs, clnsig = db2.strip().split(',')
        key2 = [chr, pos]
        value2 = [rs, geneid, clndn, clnhgvs, clnsig]
        dict_clin[','.join(key2)] = ','.join(value2)
    return dict_cos,dict_clin

def annotation_variant(dict_cos,dict_clin,variant_vcf,annotated_vcf):
    key_list = []
    key = ''
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_unmatch = 0
    var = open(variant_vcf, 'r')
    output = open(annotated_vcf, 'w')
    output.write(
        'CHR,POS,REF,ALT,QUAL,FILTER,TYPE,DP,MT,UMT,PI,THR,VMT,VMF,VSM,FORMAT,DETAIL,RS_ID,Gene_ID,CLNDN,HGVS,CLNSIG,'
        'HGNC_ID,COSMIC_ID,Mutation_Description,Feature_ID,Gene_Name,Gene_CDS_Length,Mutation_Zygosity,LOH,Mutation_Strand,'
        'HGVS.c,HGVS.p,FATHMM_Prediction,FATHMM_Score,Mutation_Somatic_Status\n')
    for line in var:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
            chrom = chrom[3:]
            detail = detail.replace(',', ';')
            # TYPE=SNP;DP=1178;MT=61;UMT=60;PI=115.34;THR=60;VMT=29;VMF=0.4833;VSM=21
            type, dp, mt, umt, pi, thr, vmt, vmf, vsm = info.split(';')
            type = type[type.find('=')+1:]
            dp = dp[dp.find('=')+1:]
            mt = mt[mt.find('=')+1:]
            umt = umt[umt.find('=')+1:]
            pi = pi[pi.find('=')+1:]
            thr = thr[thr.find('=')+1:]
            vmt = vmt[vmt.find('=')+1:]
            vmf = vmf[vmf.find('=')+1:]
            vsm = vsm[vsm.find('=')+1:]
            key = chrom + ',' + pos
            value = [chrom, pos, ref, alt, qual, filter, type, dp, mt, umt, pi, thr, vmt, vmf, vsm, format, detail]
            unmatch = 0
            # drop duplicate variant
            if key in key_list:
                continue
            if key in dict_clin:
                new = ','.join(value) + ',' + dict_clin[key] + ','
                num_in_clinvar += 1
            else:
                new = ','.join(value) + ',-,-,-,-,-,'
                unmatch += 1
            if key in dict_cos:
                new += dict_cos[key] + '\n'
                num_in_cosmic += 1
            else:
                new += '-,-,-,-,-,-,-,-,-,-,-,-,-,-\n'
                unmatch += 1
            # 2 databases both unmatch
            if unmatch == 2:
                num_unmatch += 1
                continue
            output.write(new)
        key_list.append(key)
    print ('The sample has %s variants.' % len(key_list))
    print ('%s variants in COSMIC database' % num_in_cosmic)
    print ('%s variants in Clinvar database' % num_in_clinvar)
    print ('%s variants unmatch in cosmic and clinvar.' % num_unmatch)

def main():
    time_start = time.time()
    (source, sample_name, cosmic, clinvar) = sys.argv[1:]
    variant_vcf = source + sample_name + '_variant.smCounter.cut.vcf'
    annotated_vcf = source + sample_name + '_annotated.vcf'
    dict_cos, dict_clin = read_database(cosmic,clinvar)
    annotation_variant(dict_cos,dict_clin,variant_vcf,annotated_vcf)
    print('The time use to annotate is %s minutes.' % str((time.time() - time_start) / 60))

if __name__ == '__main__':
    main()
