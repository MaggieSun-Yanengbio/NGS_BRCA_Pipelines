# INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">
# INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
# INFO=<ID=MT,Number=1,Type=Integer,Description="Total MT depth">
# INFO=<ID=UMT,Number=1,Type=Integer,Description="Filtered MT depth">
# INFO=<ID=PI,Number=1,Type=Float,Description="Variant prediction index">
# INFO=<ID=THR,Number=1,Type=Integer,Description="Variant prediction index minimum threshold">
# INFO=<ID=VMT,Number=1,Type=Integer,Description="Variant MT depth">
# INFO=<ID=VMF,Number=1,Type=Float,Description="Variant MT fraction">
# INFO=<ID=VSM,Number=1,Type=Integer,Description="Variant strong MT depth">
# FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">
# FORMAT=<ID=VF,Number=1,Type=Float,Description="Variant MT fraction, same as VMF">
# CLNDN--ClinVar`s preferred disease name for the concept specified by disease identifiers in CLNDISD.
# HGVS--Top-level(primary assembly,alt,or patch) HGVS expression.
# CLNSIG--Clinical significance for this single variant.
# Mutation_Description--Type of mutation at the amino acid level (substitution, deletion, insertion, complex, fusion, unknown etc.)
# Gene_CDS_Length--Length of the gene (base pair) counts.
# Mutation_Zygosity--Information on whether the mutation was reported to be homozygous , heterozygous or unknown within the sample.
# LOH--LOH Information on whether the gene was reported to have loss of heterozygosity in the sample: yes, no or unknown.
# Mutation_Strand--postive or negative.
# FATHMM_Prediction--Functional Analysis through Hidden Markov Models.
# FATHMM_Score--The scores are in the form of pvalues ranging from 0 to 1. Pvalues greater than 0.5 are pathogenic
# while less than 0.5 are benign. Pvalues close to 0 or 1 are the high confidence results which
# are more accurate. The results are annotated as 10 feature groups (separately for coding and
#  non coding variants) details of which can be found in the original FATHMM-MKL paper.
# Mutation_Somatic_Status--Information on whether the sample was reported to be Confirmed Somatic, Previously Reported or Variant of unknown origin.
import sys
import time
import pandas as pd

def read_database(cosmic,clinvar):
    cos = open(cosmic, 'r')
    dict_cos = {}
    cos.readline()
    for db1 in cos.readlines():
        # if not db.startswith('Gene_name'):
        Gene_name, Accession_Number, Gene_CDS_length, HGNC_ID, Sample_name, ID_sample, ID_tumour, Primary_site, \
        Site_subtype1, Site_subtype2, Site_subtype3, Primary_histology, Histology_subtype1, Histology_subtype2, \
        Histology_subtype3, Genome_wide_screen, Mutation_ID, Mutation_CDS, Mutation_AA, Mutation_Description, \
        Mutation_zygosity, LOH, GRCh, Chr, Start, End, Mutation_strand, SNP, Resistance_Mutation, FATHMM_prediction, \
        FATHMM_score, Mutation_somatic_status, Pubmed_PMID, ID_STUDY, Sample_source, Tumor_origin, Age = db1.strip().split(',')
        if Mutation_CDS is 'NS':
            continue
        if 'del' in Mutation_CDS:
            Change = Mutation_CDS[Mutation_CDS.find('del'):]
        elif 'ins' in Mutation_CDS:
            Change = Mutation_CDS[Mutation_CDS.find('ins'):]
        elif '>' in Mutation_CDS:
            Change = Mutation_CDS[Mutation_CDS.find('>')-1:]
        key1 = [Chr, Start ,Change]
        value1 = [Mutation_ID, Mutation_Description, Accession_Number, Gene_name, Gene_CDS_length,
                  Mutation_zygosity, LOH, Mutation_strand, Mutation_CDS, Mutation_AA, FATHMM_prediction, FATHMM_score,
                  Mutation_somatic_status]
        dict_cos[','.join(key1)] = ','.join(value1)
    clin = open(clinvar, 'r')
    dict_clin = {}
    clin.readline()
    for db2 in clin.readlines():
        genename, chr, start, end, geneid, alleleid, rs, pos, ref, alt, af_esp, af_exac, af_tgp, \
        clndn, clnhgvs, clnsig = db2.strip().split(',')
        if rs is 'N':
            rs = '-'
        if len(ref) == len(alt):
            change = ref + '>' + alt
        elif len(ref) > len(alt) and len(alt) == 1:
            change = 'del' + ref[1:]
        elif len(ref) < len(alt) and len(ref) == 1:
            change = 'ins' + alt[1:]
        else:
            change = 'del' + ref + 'ins' + alt
        key2 = [chr, pos, change]
        value2 = [geneid, rs, clndn, clnhgvs, clnsig]
        dict_clin[','.join(key2)] = ','.join(value2)
    return dict_cos,dict_clin


def get_info(tag):
    return tag[tag.find('=')+1:]

def annotation(dict_cos,dict_clin,variant_vcf,annotated_csv,stats_file):
    key_list = []
    key = ''
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_unmatch = 0
    var = open(variant_vcf, 'r')
    output = open(annotated_csv, 'w')
    output.write(
        'CHR,POS,REF,ALT,QUAL,FILTER,TYPE,DP,MT,UMT,PI,THR,VMT,VMF,VSM,Gene_ID,RS_ID,CLNDN,HGVS,CLNSIG,'
        'COSMIC_ID,Mutation_Description,Feature_ID,Gene_Name,Gene_CDS_Length,Mutation_Zygosity,LOH,Mutation_Strand,'
        'HGVS.c,HGVS.p,FATHMM_Prediction,FATHMM_Score,Mutation_Somatic_Status\n')
    for line in var:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
            chrom = chrom[3:]
            if len(ref) == len(alt):
                change = ref + '>' + alt
            elif len(ref) > len(alt) and len(alt) == 1:
                change = 'del' + ref[1:]
            elif len(ref) < len(alt) and len(ref) == 1:
                change = 'ins' + alt[1:]
            else:
                change = 'del' + ref + 'ins' + alt
            key = chrom + ',' + pos + ',' + change
            type, dp, mt, umt, pi, thr, vmt, vmf, vsm = map(get_info, info.split(';'))
            value = [chrom, pos, ref, alt, qual, filter, type, dp, mt, umt, pi, thr, vmt, vmf, vsm]
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
                new += '-,-,-,-,-,-,-,-,-,-,-,-,-\n'
                unmatch += 1
            # 2 databases both unmatch
            if unmatch == 2:
                num_unmatch += 1
            else:
                output.write(new)
            key_list.append(key)
    print ('The sample has %s variants.\n' % len(key_list))
    print ('%s variants in COSMIC database\n' % num_in_cosmic)
    print ('%s variants in Clinvar database\n' % num_in_clinvar)
    print ('%s variants unmatch in cosmic and clinvar\n.' % num_unmatch)
    stats_out = open(stats_file,'w')
    stats_out.write('The sample has %s variants.\n' % len(key_list))
    stats_out.write('%s variants in COSMIC database\n' % num_in_cosmic)
    stats_out.write('%s variants in Clinvar database\n' % num_in_clinvar)
    stats_out.write('%s variants unmatch in cosmic and clinvar\n.' % num_unmatch)
    stats_out.close()

def fill_table(annotated_csv,non_rs,non_cos,ref_ens):
    g2n = {}
    g2t = {}
    n2i = {}
    f1 = open(ref_ens,'r')
    for line in f1.readlines():
        l1,l2,l3 = line.strip().split(',')
        g2n[l2] = l1
        g2t[l2] = l3
        n2i[l1] = l2
    df = pd.read_csv(annotated_csv)
    subframe = df[['Gene_Name','Gene_ID','Feature_ID']]
    #for name,id,transcript in subframe.iterrows():
    for num in range(0, len(subframe)):
        #subframe.iloc[i]['Gene_Name'], subframe.iloc[i]['Gene_ID'], subframe.iloc[i]['Feature_ID']
        if subframe.iloc[num]['Gene_Name'] is '-' and subframe.iloc[num]['Feature_ID'] is '-':
            subframe.iloc[num]['Gene_Name'] = g2n[subframe.iloc[num]['Gene_ID']]
            subframe.iloc[num]['Feature_ID'] = g2t[subframe.iloc[num]['Gene_ID']]
        elif subframe.iloc[num]['Gene_ID'] is '-':
            subframe.iloc[num]['Gene_ID'] = n2i[subframe.iloc[num]['Gene_Name']]
    name = subframe['Gene_Name']
    ensg = subframe['Gene_ID']
    enst = subframe['Feature_ID']
    df.drop(labels=['Gene_Name'], axis=1, inplace=True)
    df.drop(labels=['Gene_ID'], axis=1, inplace=True)
    df.drop(labels=['Feature_ID'], axis=1, inplace=True)
    df.insert(7, 'Gene_Name', name)
    df.insert(8, 'Gene_ID', ensg)
    df.insert(9, 'Feature_ID', enst)
    df.to_csv(annotated_csv, index=False, sep=',')
    df[(True^df['RS_ID'].isin(['-']))].to_csv(non_rs, index=False, sep=',')
    df[(True^df['COSMIC_ID'].isin(['-']))].to_csv(non_cos, index=False, sep=',')

def main():
    (source, sample_name, cosmic, clinvar, ref_ens) = sys.argv[1:]
    variant_vcf = source + sample_name + '_variant.smCounter.cut.vcf'
    annotated_csv = source + sample_name + '_annotated.csv'
    non_rs = source + sample_name + '_non_rs.csv'
    non_cos = source + sample_name + '_non_cos.csv'
    stats_file = source + sample_name + '_annotate_stats.txt'
    dict_cos, dict_clin = read_database(cosmic,clinvar)
    annotation(dict_cos,dict_clin,variant_vcf,annotated_csv,stats_file)
    fill_table(annotated_csv,non_rs,non_cos,ref_ens)

if __name__ == '__main__':
    main()
