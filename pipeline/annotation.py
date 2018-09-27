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
import pymysql

base_paired = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def fetch_mysql():
    clinvar_ds = {}
    cosmic_ds = {}
    g1000_ds = {}

    connect = pymysql.connect(host = 'localhost',
                              user = 'root',
                              passwd = 'admin',
                              db = 'annotation')
    cursor = connect.cursor()
    cursor.execute('use annotation;')
    cursor.execute('select * from brca_clinvar;')
    clinvar = cursor.fetchall()
    for d1 in clinvar:
        genename, chr, start, end, gene_id, pos, ref, alt, allele_id, rs, af_esp, af_exac, af_tgp, clngn, clnhgvs, \
        clnsig, so, molecular_consequence, feature_id, exon, hgvs_c, hgvs_p = d1
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
        k1 = [chr, pos, change]
        v1 = [gene_id, rs, clngn, clnhgvs, clnsig, so, molecular_consequence]
        clinvar_ds[','.join(k1)] = ','.join(v1)

    cursor.execute('select * from cosmic;')
    cosmic = cursor.fetchall()
    for d2 in cosmic:
        genename, accessing_number, cds_length, hgnc_id, sample_name, sample_id, tumor_id, primary_site, site_sub1, \
        site_sub2, site_sub3, primary_histology, hist_sub1, hist_sub2, hist_sub3, genome_wide_screen, mutation_id, \
        mutation_cds, mutation_aa, mutation_description, mutation_zygosity, loh, grch, chr, start, end, hgvs, \
        mutation_strand, snp, resistance, prediction, score, mutation_somatic_status, pubmed, id_study, sample_type, \
        tumor_origin, age = d2
        if mutation_cds is 'NS':
            continue
        if 'del' in mutation_cds:
            change = mutation_cds[mutation_cds.find('del'):]
        elif 'ins' in mutation_cds:
            change = mutation_cds[mutation_cds.find('ins'):]
        elif '>' in mutation_cds:
            change = mutation_cds[mutation_cds.find('>')-1:]
        k2 = [chr, start ,change]
        v2 = [mutation_id, mutation_description, accessing_number, genename, mutation_zygosity, loh, prediction,
                  score, mutation_somatic_status, pubmed]
        cosmic_ds[','.join(k2)] = ','.join(v2)

    cursor.execute('select * from g1000;')
    g1000 = cursor.fetchall()
    for d3 in g1000:
        genename, chr, start, end, variant_type, ref, alt, rs, eas_af, eur_af, amr_af, sas_af, afr_af, hgvs = d3
        if variant_type == 'SNV':
            change = ref + '>' + alt
        elif variant_type == 'insertion':
            change = 'ins' + alt
        elif variant_type == 'deletion':
            change = 'del' + ref
        else:
            pass
        k3 = [chr, start, change]
        v3 = [genename, rs, eas_af, eur_af, amr_af, sas_af, afr_af]
        g1000_ds[','.join(k3)] = ','.join(v3)

    cursor.close()
    connect.commit()
    connect.close()

    return cosmic_ds, clinvar_ds, g1000_ds


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
        'CHR,POS,REF,ALT,FILTER,TYPE,DP,MT,UMT,VMT,VMF,TLOD,NLOD,Gene_ID,RS_ID,CLNDN,HGVS,CLNSIG,SO,Molecular_Consequence_Clinvar,'
        'COSMIC_ID,Mutation_Consequence_Cosmic,Feature_ID,Gene_Name,Mutation_Zygosity,LOH,'
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