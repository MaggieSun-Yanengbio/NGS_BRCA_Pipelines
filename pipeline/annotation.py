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

base_paired = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def read_database(cosmic,clinvar,g1000):
    cos = open(cosmic, 'r')
    dict_cos = {}
    cos.readline()
    for db1 in cos.readlines():
        # if not db.startswith('Gene_name'):
        gene_name, accession_number, gene_cds_length, hgnc_id, sample_name, id_sample, id_tumour, primary_site, \
        site_subtype1, site_subtype2, site_subtype3, primary_histology, histology_subtype1, histology_subtype2, \
        histology_subtype3, genome_wide_screen, mutation_id, mutation_cds, mutation_aa, mutation_description, \
        mutation_zygosity, loh, grch, chr, start, end, mutation_strand, snp, resistance_mutation, fathmm_prediction, \
        fathmm_score, mutation_somatic_status, pubmed_pmid, id_study, sample_source, tumor_origin, age = db1.strip().split(',')
        if mutation_cds is 'NS':
            continue
        if 'del' in mutation_cds:
            change = mutation_cds[mutation_cds.find('del'):]
        elif 'ins' in mutation_cds:
            change = mutation_cds[mutation_cds.find('ins'):]
        elif '>' in mutation_cds:
            change = mutation_cds[mutation_cds.find('>')-1:]
        key1 = [chr, start ,change]
        value1 = [mutation_id, mutation_description, accession_number, gene_name,
                  mutation_zygosity, loh, fathmm_prediction, fathmm_score,
                  mutation_somatic_status, pubmed_pmid]
        dict_cos[','.join(key1)] = ','.join(value1)

    clin = open(clinvar, 'r')
    dict_clin = {}
    clin.readline()
    for db2 in clin.readlines():
        genename, chr, start, end, geneid, pos, ref, alt, alleleid, rs, af_esp, af_exac, af_tgp, \
        clndn, clnhgvs, clnsig, so, molecular_consequence = db2.strip().split(',')
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
        value2 = [geneid, rs, clndn, clnhgvs, clnsig, so, molecular_consequence]
        dict_clin[','.join(key2)] = ','.join(value2)

    genomes1000 = open(g1000, 'r')
    dict_g1000 = {}
    genomes1000.readline()

    for db3 in genomes1000.readlines():
        genename1, chr1, start1, end1, variant_type, ref1, alt1, rs1, eas_af, eur_af, amr_af, sas_af, afr_af, hgvs1 = db3.strip().split(',')
        if variant_type == 'SNV':
            change1 = ref1 + '>' + alt1
        elif variant_type == 'insertion':
            change1 = 'ins' + alt1
        elif variant_type == 'deletion':
            change1 = 'del' + ref1
        else:
            pass
        key3 = [chr1, start1, change1]
        value3 = [genename1, rs1, eas_af, eur_af, amr_af, sas_af, afr_af]
        dict_g1000[','.join(key3)] = ','.join(value3)
    return dict_cos,dict_clin,dict_g1000

def get_info(tag):
    return tag[tag.find('=')+1:]


def define_hgvs(chr, pos, ref, alt):
    chr_to_version = {'1': '10', '2': '11', '3': '11', '4': '11', '5': '9', '6': '11', '7': '13', '8': '10', '9': '11',
                      '10': '10', '11': '9', '12': '11', '13': '10', '14': '8', '15': '9', '16': '9', '17': '10',
                      '18': '9', '19': '9', '20': '10', '21': '8', '22': '10', 'X': '10'}
    version = chr_to_version[chr]
    if chr != 'X' and int(chr) < 10:
        chr = '0' + chr
    if len(ref) == 1 and len(alt) == 1:
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + pos + ref + '>' + alt
    if len(ref) > 1 and len(alt) == 1:
        deletion = ref[1:]
        if len(deletion) == 1:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+1) + 'del' + deletion
        else:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+1) + '_' + str(int(pos)+len(deletion)) + 'del' + deletion
    if len(ref) == 1 and len(alt) > 1:
        insertion = alt[1:]
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)) + '_' + str(int(pos)+1) + 'ins' + insertion

    return hgvs


def annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,stats_file):
    key_list = []
    key = ''
    change1 = ''
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_in_g1000 = 0
    num_unmatch = 0
    var = open(variant_vcf, 'r')
    output = open(annotated_csv, 'w')
    output.write(
        'CHR,POS,REF,ALT,FILTER,TYPE,DP,MT,UMT,VMT,VMF,TLOD,NLOD,Gene_ID,RS_ID,CLNDN,HGVS,CLNSIG,SO,Molecular_Consequence,'
        'COSMIC_ID,Mutation_Description,Feature_ID,Gene_Name,Mutation_Zygosity,LOH,'
        'FATHMM_Prediction,FATHMM_Score,Mutation_Somatic_Status,PMID,Gene_Name1,RS_ID1,EAS_AF,EUR_AF,AMR_AF,'
        'SAS_AF,AFR_AF\n')
    for line in var:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
            if filter == 'PASS':
                chrom = chrom[3:]
                if len(ref) == len(alt):
                    change = ref + '>' + alt
                    change1 = base_paired[ref] + '>' + base_paired[alt]
                elif len(ref) > len(alt) and len(alt) == 1:
                    change = 'del' + ref[1:]
                elif len(ref) < len(alt) and len(ref) == 1:
                    change = 'ins' + alt[1:]
                else:
                    change = 'del' + ref + 'ins' + alt
                key = chrom + ',' + pos + ',' + change
                key1 = chrom + ',' + pos + ',' + change1
                type, dp, mt, umt, pi, thr, vmt, vmf, vsm, tlod, nlod = map(get_info, info.split(';'))
                value = [chrom, pos, ref, alt, filter, type, dp, mt, umt, vmt, vmf, tlod, nlod]
                unmatch = 0
                # drop duplicate variant
                if key in key_list:
                    continue
                #if filter == 'PASS':
                if key in dict_clin:
                    new = ','.join(value) + ',' + dict_clin[key] + ','
                    num_in_clinvar += 1
                else:
                    new = ','.join(value) + ',-,-,-,-,-,' + define_hgvs(chrom, pos, ref, alt) + ',-,'
                    unmatch += 1
                if key in dict_cos:
                    new += dict_cos[key] + ','
                    num_in_cosmic += 1
                elif key1 in dict_cos:
                    new += dict_cos[key1] + ','
                    num_in_cosmic += 1
                else:
                    new += '-,-,-,-,-,-,-,-,-,-,'
                    unmatch += 1
                if key in dict_g1000:
                    new += dict_g1000[key] + '\n'
                    num_in_g1000 += 1
                else:
                    new += '-,-,-,-,-,-,-\n'
                    unmatch += 1
                # 3 databases both unmatch
                if unmatch == 3:
                    num_unmatch += 1
                else:
                    output.write(new)
                key_list.append(key)
    print ('The sample has %s variants.' % len(key_list))
    print ('%s variants in cosmic database.' % num_in_cosmic)
    print ('%s variants in clinvar database.' % num_in_clinvar)
    print ('%s variants in 1000genomes database.' % num_in_g1000)
    print ('%s variants unmatch in cosmic,clinvar and 1000genomes.\n' % num_unmatch)
    stats_out = open(stats_file,'w')
    stats_out.write('The sample has %s variants.' % len(key_list))
    stats_out.write('%s variants in cosmic database.' % num_in_cosmic)
    stats_out.write('%s variants in clinvar database.' % num_in_clinvar)
    stats_out.write('%s variants in 1000genomes database.' % num_in_g1000)
    stats_out.write('%s variants unmatch in cosmic,clinvar and 1000genomes.\n' % num_unmatch)
    stats_out.close()

def fill_table(annotated_csv, ref_ens):
    n2g = {}
    g2n = {}
    for line in open(ref_ens, 'r').readlines():
        name, ensg = line.strip().split(',')
        n2g[name] = ensg
        g2n[ensg] = name

    df = pd.read_csv(annotated_csv)
    subframe = df[['Gene_Name','Gene_ID', 'Gene_Name1','RS_ID','RS_ID1']]
    for num in range(0, len(subframe)):
        if subframe.iloc[num, 0] is '-' and subframe.iloc[num, 1] is '-':
            name = subframe.iloc[num, 2]
            ensg = n2g[name]
            subframe.iloc[num, 0] = name
            subframe.iloc[num, 1] = ensg
        elif subframe.iloc[num, 0] is '-' and subframe.iloc[num, 2] is '-':
            ensg = subframe.iloc[num, 1]
            name = g2n[ensg]
            subframe.iloc[num, 0] = name
        elif subframe.iloc[num, 1] is '-' and subframe.iloc[num, 2] is '-':
            name = subframe.iloc[num, 0]
            ensg = n2g[name]
            subframe.iloc[num, 1] = ensg
        elif subframe.iloc[num, 0] is '-':
            name = subframe.iloc[num, 2]
            subframe.iloc[num, 0] = name
        elif subframe.iloc[num, 1] is '-':
            name = subframe.iloc[num, 0]
            ensg = n2g[name]
            subframe.iloc[num, 1] = ensg
        if subframe.iloc[num, 3] is '-' and subframe.iloc[num, 4] is not '-':
            subframe.iloc[num, 3] = subframe.iloc[num, 4]
    name = subframe['Gene_Name']
    ensg = subframe['Gene_ID']
    rs = subframe['RS_ID']
    df.drop(labels=['Gene_Name'], axis=1, inplace=True)
    df.drop(labels=['Gene_ID'], axis=1, inplace=True)
    df.drop(labels=['Gene_Name1'], axis=1, inplace=True)
    df.drop(labels=['RS_ID'], axis=1, inplace=True)
    df.drop(labels=['RS_ID1'], axis=1, inplace=True)
    df.insert(7, 'Gene_Name', name)
    df.insert(8, 'Gene_ID', ensg)
    df.insert(7, 'RS_ID', rs)
    df.to_csv(annotated_csv, index=False, sep=',')