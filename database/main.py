import sys
import os
import re
import pymysql
from functools import reduce
from subprocess import Popen, PIPE
from collections import defaultdict


base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
genes = ['BRCA1', 'BRCA2']
# data from genecards
genomic_locations = {'BRCA1':['17', '41196312', '41277500'], 'BRCA2':['13', '32889611', '32973809']}
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
               '18', '19', '20', '21', '22', 'X']

def define_hgvs(chr, pos, ref, alt):
    chr_to_version = {'1': '10', '2': '11', '3': '11', '4': '11', '5': '9', '6': '11', '7': '13', '8': '10', '9': '11',
                      '10': '10', '11': '9', '12': '11', '13': '10', '14': '8', '15': '9', '16': '9', '17': '10',
                      '18': '9', '19': '9', '20': '10', '21': '8', '22': '10', '23': '10'}
    # chr = chr[3:]
    ref = ref.upper()
    alt = alt.upper()
    if chr == 'X':
        chr = '23'
    version = chr_to_version[chr]
    if int(chr) < 10:
        chr = '0' + chr
    if len(ref) == 1 and len(alt) == 1:
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + pos + ref + '>' + alt
    elif len(ref) > 1 and len(alt) == 1:
        deletion = ref[1:]
        if len(deletion) == 1:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos) + 1) + 'del' + deletion
        else:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos) + 1) + '_' + str(
                int(pos) + len(deletion)) + 'del' + deletion
    elif len(ref) == 1 and len(alt) > 1:
        insertion = alt[1:]
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + pos + '_' + str(int(pos) + 1) + 'ins' + insertion
    else:
        deletion = ref
        insertion = alt
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + pos + '_' + str(int(pos) + len(deletion) - 1) + 'del'\
               + deletion + 'ins' + insertion
    return hgvs


def refbase(samtools, genome, chr, position):
    stdout, stderr = Popen([samtools, 'faidx', genome, 'chr{}:{}-{}'.format(chr, position, position)],
                           stdout=PIPE).communicate()
    base = str(stdout, encoding='utf-8').split('\n')[1]
    return base


def reverse_complementary(seq):
    seq = [base_pair[i] for i in seq]
    return reduce(lambda x, y: y + x, seq)


def login_mysql(host, user, passwd, database):
    conn = pymysql.connect(
        host = host,
        user = user,
        passwd = passwd,
        db = database
    )
    cur = conn.cursor()
    return conn, cur


def check_table(cursor, list):
    # check row
    sql = "select * from mutations where HGVS='{}';".format(list[0][1])
    cursor.execute(sql)
    result1 = cursor.fetchone()
    if result1:
        # check column
        for i in list[6:]:
            sql = "select * from mutations where HGVS='{}' and {}='{}'".format(list[0][1], i[0], i[1])
            cursor.execute(sql)
            result2 = cursor.fetchone()
            if result2:
                pass
            else:
                sql = "update mutations set {}='{}' where HGVS='{}'".format(i[0], i[1], list[0][1])
                # sql = "insert into mutations({}) values('{}');".format(i[0], i[1])
                cursor.execute(sql)
    else:
        sql = "insert into mutations(HGVS) values('{}');".format(list[0][1])
        cursor.execute(sql)
        for i in list[1:]:
            sql = "update mutations set {}='{}' where HGVS='{}'".format(i[0], i[1], list[0][1])
            # sql = "insert into mutations({}) values('{}');".format(i[0], i[1])
            cursor.execute(sql)


def brca1_cdna_to_genomes(exons_position, gene, exon, c_position):
    if int(exon) > 3:
        exon = str(int(exon) - 1)
    exon_start, exon_end, pass_length = exons_position[gene][exon]
    if exon == '2':
        true_pos = str(41277500 - 1368 - 19 - (int(c_position) - 1))
    else:
        true_pos = str(41277500 - (int(exon_start) - 1) - ((int(c_position) - 1) - (int(pass_length) - (213 + 19))))
    return true_pos


def brca2_cdna_to_genomes(exons_position, gene, exon, c_position):

    exon_start, exon_end, pass_length = exons_position[gene][exon]
    if exon == '2':
        true_pos = str(int(exon_start) + 39 + (int(c_position) - 1))
    else:
        true_pos = str(int(exon_start) + ((int(c_position) - 1) - (int(pass_length) - (188 + 39))))
    return true_pos


def update_clinvar(vcf, cursor):
    with open(vcf, 'r') as filein:
        for row in filein:
            tag_to_content = defaultdict(str)
            if not row.startswith('#'):
                chr, pos, id, ref, alt, qual, filter, info = row.strip().split('\t')
                if 'GENEINFO' in info and chr in chromosomes and len(ref) < 500 and len(alt) < 500:
                    for tag in info.strip().split(';'):
                        tag_to_content[tag[:tag.rfind('=')]] = tag[tag.rfind('=')+1:]
                    if tag_to_content['CLNDN'] == '':
                        tag_to_content['CLNDN'] = '-'
                    else:
                        tag_to_content['CLNDN'] = tag_to_content['CLNDN'].replace(',', ';')
                    if tag_to_content['CLNSIG'] == '':
                        tag_to_content['CLNSIG'] = '-'
                    else:
                        tag_to_content['CLNSIG'] = tag_to_content['CLNSIG'].replace(',', ';')
                    if 'MC' in tag_to_content:
                        so, molecular_consequence = tag_to_content['MC'].split(',')[0].split('|')
                        so = so[3:]
                    else:
                        so, molecular_consequence = ['-', '-']

                    geneinfo = tag_to_content['GENEINFO']
                    if '|' not in tag_to_content['GENEINFO']:
                        genename = geneinfo[geneinfo.rfind('=')+1:geneinfo.rfind(':')]
                        hgvs = define_hgvs(chr, pos, ref, alt)
                        if genename in genes:
                            check_table(cursor, [['HGVS', hgvs], ['GENE', genename], ['CHR', chr], ['POSITION', pos],
                                                 ['REF', ref], ['ALT', alt], ['CLNDN_CLINVAR', tag_to_content['CLNDN']],
                                                 ['SIGNIFICANCE_CLINVAR', tag_to_content['CLNSIG']], ['SO_CLINVAR', so],
                                                 ['MUTATION_DESCRIPTION_CLINVAR', molecular_consequence]])
                    else:
                        for gene in geneinfo[geneinfo.find('=')+1:].split('|'):
                            genename = gene[:gene.rfind(':')]
                            hgvs = define_hgvs(chr, pos, ref, alt)
                            if genename in genes:
                                check_table(cursor, [['HGVS', hgvs], ['GENE', genename], ['CHR', chr], ['POSITION', pos],
                                                     ['REF', ref], ['ALT', alt], ['CLNDN_CLINVAR', tag_to_content['CLNDN']],
                                                     ['SIGNIFICANCE_CLINVAR', tag_to_content['CLNSIG']],
                                                     ['SO_CLINVAR', so],
                                                     ['MUTATION_DESCRIPTION_CLINVAR', molecular_consequence]])
    return 0


def update_cosmic(vcf, cursor):
    hgvs_set = []
    with open(vcf, 'r') as filein:
        filein.readline()
        for row in filein:
            genename, accessing_number, cds_length, hgnc_id, sample_name, sample_id, tumor_id, primary_site, site_sub1, \
            site_sub2, site_sub3, primary_histology, hist_sub1, hist_sub2, hist_sub3, genome_wide_screen, mutation_id, \
            mutation_cds, mutation_aa, mutation_description, mutation_zygosity, log, grch, pos, mutation_strand, snp, \
            resistance, prediction, score, mutation_somatic_status, pubmed, id_study, sample_type, tumor_origin, age \
                = row.split('\t')
            if 'ENST' in genename:
                genename = genename[:genename.find('_')]
            if genename in genes and mutation_cds != '?' and mutation_aa != '?' and mutation_description != 'Unknown' and pos.strip() != '':
                chr = pos[:pos.find(':')]
                start = pos[pos.find(':') + 1:pos.find('-')]
                end = pos[pos.find('-') + 1:]
                if prediction.strip() == '':
                    prediction = '-'
                if score.strip() == '':
                    score = '-'
                if pubmed.strip() == '':
                    pubmed = '-'
                change = re.split('[\d+]', mutation_cds)[-1]
                # init
                ref = False
                alt = False
                if '>' in mutation_cds:
                    ref, alt = change.split('>')
                    position = start
                    if len(ref) >= 1 and len(alt) >= 1:
                        if genename == 'BRCA1':
                            ref, alt = reverse_complementary(ref), reverse_complementary(alt)
                    else:
                        ref = False
                        alt = False
                elif 'del' in mutation_cds:
                    deletion = mutation_cds[mutation_cds.find('del') + 3:]
                    position = str(int(start) - 1)
                    if deletion.isalpha():
                        last_base = refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chr, position)
                        if genename == 'BRCA1':
                            ref = last_base + reverse_complementary(deletion)
                        else:
                            ref = last_base + deletion
                        alt = last_base
                elif 'ins' in mutation_cds:
                    insertion = mutation_cds[mutation_cds.find('ins') + 3:]
                    position = start
                    if insertion.isalpha():
                        last_base = refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chr, position)
                        ref = last_base
                        if genename == 'BRCA1':
                            alt = last_base + reverse_complementary(insertion)
                        else:
                            alt = last_base + insertion
                if ref and alt:
                    hgvs = define_hgvs(chr, position, ref, alt)
                    if hgvs and hgvs not in hgvs_set:
                        check_table(cursor, [['HGVS', hgvs], ['GENE', genename], ['CHR', chr], ['POSITION', position],
                                             ['REF', ref], ['ALT', alt], ['COSMID_COSMIC', mutation_id],
                                             ['FEATURE_ID_COSMIC', accessing_number], ['CDS_LENGTH_COSMIC', cds_length],
                                             ['HGVS_C_COSMIC', mutation_cds], ['HGVS_P_COSMIC', mutation_aa],
                                             ['MUTATION_DESCRIPTION_COSMIC', mutation_description],
                                             ['FATHMM_PREDICT_COSMIC', prediction],
                                             ['FATHMM_SCORE_COSMIC', score], ['PMID_COSMIC', pubmed]])
                        hgvs_set.append(hgvs)
    return 0


def update_1000genomes(vcf, cursor):
    brca1 = genomic_locations['BRCA1']
    brca2 = genomic_locations['BRCA2']
    with open(vcf, 'r') as filein:
        for row in filein:
            tag_to_content = defaultdict(str)
            if not row.startswith('#'):
                chrom, pos, rs, ref, alt, qual, flt, info = row.split('\t')[:8]
                for tag in info.strip().split(';'):
                    tag_to_content[tag[:tag.rfind('=')]] = tag[tag.rfind('=') + 1:]
                hgvs = define_hgvs(chrom, pos, ref, alt)
                if chrom == brca1[0] and int(pos) >= int(brca1[1]) and int(pos) <= int(brca1[2]):
                    check_table(cursor, [['HGVS', hgvs], ['GENE', 'BRCA1'], ['CHR', chr], ['POSITION', pos],
                                         ['REF', ref], ['ALT', alt], ['RS_1000GENOMES', rs],
                                         ['AF_1000GENOMES', tag_to_content['AF']],
                                         ['EAS_AF_1000GENOMES', tag_to_content['EAS_AF']],
                                         ['AMR_AF_1000GENOMES', tag_to_content['AMR_AF']],
                                         ['AFR_AF_1000GENOMES', tag_to_content['AFR_AF']],
                                         ['EUR_AF_1000GENOMES', tag_to_content['EUR_AF']],
                                         ['SAS_AF_1000GENOMES', tag_to_content['SAS_AF']], ])
                elif chrom == brca2[0] and int(pos) >= int(brca2[1]) and int(pos) <= int(brca2[2]):
                    check_table(cursor, [['HGVS', hgvs], ['GENE', 'BRCA2'], ['CHR', chr], ['POSITION', pos],
                                         ['REF', ref], ['ALT', alt], ['RS_1000GENOMES', rs],
                                         ['AF_1000GENOMES', tag_to_content['AF']],
                                         ['EAS_AF_1000GENOMES', tag_to_content['EAS_AF']],
                                         ['AMR_AF_1000GENOMES', tag_to_content['AMR_AF']],
                                         ['AFR_AF_1000GENOMES', tag_to_content['AFR_AF']],
                                         ['EUR_AF_1000GENOMES', tag_to_content['EUR_AF']],
                                         ['SAS_AF_1000GENOMES', tag_to_content['SAS_AF']], ])

    return 0


def update_dbsnp(vcf, cursor):
    with open(vcf, 'r') as filein:
        for row in filein:
            if not row.startswith('#') and 'GENEINFO' in row:
                genelist = []
                chrom, pos, rs, ref, alt, qual, flt, info = row.strip().split('\t')
                tag_to_content = {i[:i.find('=')]: i[i.find('=') + 1:] for i in info.split(';') if '=' in i}
                for g in tag_to_content['GENEINFO'].split('|'):
                    genelist.append(g.split(':')[0])
                for gene in genes:
                    if gene in genelist:
                        for sub_alt in alt.split(','):
                            hgvs = define_hgvs(chrom, pos, ref, sub_alt)
                            check_table(cursor, [['HGVS', hgvs], ['GENE', gene], ['CHR', chrom], ['POSITION', pos],
                                                ['REF', ref], ['ALT', sub_alt], ['RS_DBSNP', rs]])
    return 0


def update_evs(vcf, cursor):
    with open(vcf, 'r') as filein:
        for row in filein:
            if not row.startswith('#'):
                chrom, pos, id, ref, alt, qual, flt, info = row.strip().split('\t')
                dbsnp, ea_ac, aa_ac, tac, maf, gts, ea_gtc, aa_gtc, gtc, dp, gl, cp, cg, aa, ca, exome_chip, gwas_pubmed, \
                fg, hgvs_c, hgvs_protein, cds_sizes, gs, ph, ea_age, aa_age, grch38 = info.split(';')
                for gene in genes:
                    if gene in gl:
                        hgvs = define_hgvs(chrom, pos, ref, alt)
                        check_table(cursor, [['HGVS', hgvs], ['GENE', gene], ['CHR', chrom], ['POSITION', pos],
                                            ['REF', ref], ['ALT', alt], ['MUTATION_DESCRIPTION_EVS', fg],
                                            ['HGVS_C_EVS', hgvs_c], ['HGVS_P_EVS', hgvs_protein], ['CDS_LENGTH_EVS', cds_sizes],
                                            ['POLYPHEN2_PREDICT_EVS', ph], ['GRCH38_POSITION_EVS', grch38]])
    return 0


def update_exac(vcf, cursor):
    with open(vcf, 'r') as filein:
        for row in filein:
            if not row.startswith('#'):
                chrom, pos, id, ref, alt, qual, flt, info = row.strip().split('\t')
                if flt == 'PASS':
                    tag_to_content = {i[:i.find('=')]: i[i.find('=') + 1:] for i in info.split(';') if '=' in i}
                    for gene in genes:
                        if gene in tag_to_content['CSQ']:
                            for a in alt.split(','):
                                hgvs = define_hgvs(chrom, pos, ref, a)
                                check_table(cursor, [['HGVS', hgvs], ['GENE', gene], ['CHR', chrom], ['POSITION', pos],
                                                     ['REF', ref], ['ALT', a], ['AF_EXAC', tag_to_content['AF']],
                                                     ['AN_EAS_EXAC', tag_to_content['AN_EAS']]])
    return 0


def update_utahdb_brca(csv, exons_position, cursor):
    brca_exons = defaultdict(dict)
    for row in open(exons_position, 'r'):
        if not row.startswith('#'):
            gene, exon, start, end, length = row.strip().split(',')
            brca_exons[gene][exon] = [start, end, length]
    for row in open(csv, 'r'):
        if not row.startswith('#'):
            # print (row)
            gene, exon, tp, change, pro, classfication, post, ref, ref2 = row.strip().split('\t')
            if gene == 'BRCA2' and 'Exon' in exon:
                chrom = '13'
                exon = exon.split(' ')[1]
                if '+' in change or '-' in change:
                    continue
                if '>' in change:
                    p = change[2:change.find('>') - 1]
                    true_pos = brca2_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('>') - 1:]
                    ref, alt = c.split('>')

                elif 'del' in change and 'ins' in change:
                    p = change[2:change.find('del')].split('_')[0]
                    true_pos = brca2_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('del'):]
                    ref, alt = c.split('ins')
                    ref = ref[3:]

                elif 'del' in change:
                    p = change[2:change.find('del')].split('_')[0]
                    true_pos = brca2_cdna_to_genomes(brca_exons, gene, exon, int(p)-1)
                    c = change[change.find('del'):]
                    last_base = refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chrom, true_pos)
                    ref = last_base + c[3:]
                    alt = last_base

                elif 'dup' in change:
                    p = change[2:change.find('dup')].split('_')[-1]
                    true_pos = brca2_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('dup'):]
                    ref = c[-1]
                    alt = ref + c[3:]

                elif 'ins' in change:
                    p = change[2:change.find('ins')].split('_')[0]
                    true_pos = brca2_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('ins'):]
                    last_base = refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chrom, true_pos)
                    ref = last_base
                    alt = last_base + c[3:]
                if ref.isalpha() and alt.isalpha():
                    # print (chrom, true_pos, ref, alt)
                    hgvs = define_hgvs(chrom, true_pos, ref, alt)
                    check_table(cursor, [['HGVS', hgvs], ['GENE', gene], ['CHR', chrom], ['POSITION', true_pos],
                                                         ['REF', ref], ['ALT', alt], ['SIGNIFICANCE_UTAHDB', classfication],
                                                         ['REFERENCE1_UTAHDB', ref], ['REFERENCE2_UTAHDB', ref2]])

            elif gene == 'BRCA1' and 'Exon' in exon:
                chrom = '17'
                exon = exon.split(' ')[1]
                if '+' in change or '-' in change:
                    continue
                if '>' in change:
                    p = change[2:change.find('>') - 1]
                    true_pos = brca1_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('>') - 1:]
                    ref, alt = c.split('>')
                    # ref, alt = [base_pair[i] for i in c.split('>')]

                elif 'del' in change and 'ins' in change:
                    p = change[2:change.find('del')].split('_')[-1]
                    true_pos = brca1_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('del'):]
                    ref, alt = c.split('ins')
                    ref = ref[3:]
                    # ref = reverse_complementary(ref[3:])
                    # alt = reverse_complementary(alt)

                elif 'del' in change:
                    p = change[2:change.find('del')].split('_')[-1]
                    true_pos = brca1_cdna_to_genomes(brca_exons, gene, exon, int(p) + 1)
                    c = change[change.find('del'):]
                    last_base = base_pair[refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chrom, true_pos)]
                    # print (true_pos, last_base,refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chrom, true_pos))
                    # ref = last_base + reverse_complementary(c[3:])
                    ref = c[3:] + last_base
                    alt = last_base

                elif 'dup' in change:
                    p = change[2:change.find('dup')].split('_')[0]
                    true_pos = brca1_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('dup'):]
                    ref = c[3]
                    alt = c[3:] + ref

                elif 'ins' in change:
                    p = change[2:change.find('ins')].split('_')[-1]
                    true_pos = brca1_cdna_to_genomes(brca_exons, gene, exon, p)
                    c = change[change.find('ins'):]
                    last_base = base_pair[refbase('samtools', '/home/administrator/source/ucsc.hg19.fasta', chrom, true_pos)]
                    ref = last_base
                    alt = c[3:] + last_base

                if ref.isalpha() and alt.isalpha():
                    # print (ref, alt)
                    ref = reverse_complementary(ref)
                    alt = reverse_complementary(alt)
                    # print (chrom, true_pos, ref, alt)
                    hgvs = define_hgvs(chrom, true_pos, ref, alt)
                    check_table(cursor, [['HGVS', hgvs], ['GENE', gene], ['CHR', chrom], ['POSITION', true_pos],
                                                         ['REF', ref], ['ALT', alt], ['SIGNIFICANCE_UTAHDB', classfication],
                                                         ['REFERENCE1_UTAHDB', ref], ['REFERENCE2_UTAHDB', ref2]])
    return 0


def update_cngb(csv, cursor):
    with open(csv, 'r') as filein:
        for row in filein:
            if not row.startswith('#'):
                chrom, pos, ref, alt, gene, rs, hpvs_p, annotation, clnsig, alt_num, num, vf = row.strip().split(',')
                if gene in genes:
                    hgvs = define_hgvs(chrom, pos, ref, alt)
                    check_table(cursor, [['HGVS', hgvs], ['GENE', gene], ['CHR', chrom], ['POSITION', pos],
                                         ['REF', ref], ['ALT', alt], ['AF_CNGB', vf]])
    return 0


def main():

    connect, cursor = login_mysql('localhost', 'root', 'admin', 'annotation')
    # # BRCA1 and BRCA2 exons position table
    # exons = '/home/administrator/database/mutation_database/exons_position.csv'
    # # input: database + vcf path
    # parameters = sys.argv[1:]
    # modes = [parameters[i:i + 2] for i in range(0, len(parameters), 2)]
    # for mode in modes:
    #     if mode[0] == 'clinvar':
    #         update_clinvar(mode[1], cursor)
    #     elif mode[0] == 'cosmic':
    #         update_cosmic(mode[1], cursor)
    #     elif mode[0] == '1000genomes':
    #         update_1000genomes(mode[1], cursor)
    #     elif mode[0] == 'dbsnp':
    #         update_dbsnp(mode[1], cursor)
    #     elif mode[0] == 'evs':
    #         update_evs(mode[1], cursor)
    #     elif mode[0] == 'exac':
    #         update_exac(mode[1], cursor)
    #     elif mode[0] == 'utahdb':
    #         update_utahdb_brca(mode[1], exons, cursor)
    #     elif mode[0] == 'cngb':
    #         update_cngb(mode[1], cursor)
    #     else:
    #         print ('Error input.')

    vcf = '/home/administrator/database/mutation_database/clinvar/clinvar_20180729.vcf'
    update_clinvar(vcf, cursor)
    vcf = '/home/administrator/database/mutation_database/cosmic/CosmicMutantExport_20180820.tsv'
    update_cosmic(vcf, cursor)
    vcf = '/home/administrator/database/mutation_database/1000genomes/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'
    update_1000genomes(vcf, cursor)
    vcf = '/home/administrator/database/mutation_database/1000genomes/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'
    update_1000genomes(vcf, cursor)
    # vcf = '/home/administrator/database/mutation_database/dbsnp/All_20170710.vcf'
    # update_dbsnp(vcf, cursor)
    vcf = '/home/administrator/database/mutation_database/evs/ESP6500SI-V2-SSA137.GRCh38-liftover.chr13.snps_indels.vcf'
    update_evs(vcf, cursor)
    vcf = '/home/administrator/database/mutation_database/evs/ESP6500SI-V2-SSA137.GRCh38-liftover.chr17.snps_indels.vcf'
    update_evs(vcf, cursor)
    vcf = '/home/administrator/database/mutation_database/exac/exac_brca.vcf'
    update_exac(vcf, cursor)
    csv = '/home/administrator/database/mutation_database/arup/arup_brca.tsv'
    exons = '/home/administrator/database/mutation_database/exons_position.csv'
    update_utahdb_brca(csv, exons, cursor)
    csv = '/home/administrator/database/mutation_database/cngb/BRCA.csv'
    update_cngb(csv, cursor)
    cursor.close()
    connect.commit()
    connect.close()


if __name__ == '__main__':
    main()