# clinvar数据整理
from collections import defaultdict

clinvar_table = open('/home/administrator/database/mutation_database/clinvar/clinvar_table.csv', 'w')
clinvar_non_genename = open('/home/administrator/database/mutation_database/clinvar/clinvar_non_genename.vcf', 'w')

for row in open('/home/administrator/database/mutation_database/clinvar/clinvar_20180729.vcf','r'):
    tag_to_content = defaultdict(str)
    if not row.startswith('#'):
        chr, pos, id, ref, alt, qual, filter, info = row.strip().split('\t')
        if 'GENEINFO' in info:
            for tag in info.strip().split(';'):
                tag_to_content[tag[:tag.rfind('=')]] = tag[tag.rfind('=')+1:]
            tag_to_content['CLNDN'] = tag_to_content['CLNDN'].replace(',',';')
            tag_to_content['CLNSIG'] = tag_to_content['CLNSIG'].replace(',',';')
            if 'MC' in tag_to_content:
                so, molecular_consequence = tag_to_content['MC'].split(',')[0].split('|')
            else:
                so, molecular_consequence = ['N', 'N']

            geneinfo = tag_to_content['GENEINFO']
            if '|' not in tag_to_content['GENEINFO']:
                genename = geneinfo[geneinfo.rfind('=')+1:geneinfo.rfind(':')]
                join_list = [genename, pos, ref, alt , tag_to_content['ALLELEID'], tag_to_content['RS'],
                                    tag_to_content['AF_ESP'], tag_to_content['AF_EXAC'], tag_to_content['AF_TGP'],
                                    tag_to_content['CLNDN'], tag_to_content['CLNHGVS'], tag_to_content['CLNSIG'],
                                    so[3:], molecular_consequence]
                idx = 0
                for ele in join_list:
                    if ele == '':
                        join_list[idx] = 'N'
                    idx += 1
                new_row = ','.join(join_list)+'\n'
                clinvar_table.write(new_row)
            else:
                for gene in geneinfo.split('|'):
                    genename = gene[:gene.rfind(':')]
                    join_list = [genename, chr, pos, ref, alt, tag_to_content['ALLELEID'], tag_to_content['RS'],
                                        tag_to_content['AF_ESP'], tag_to_content['AF_EXAC'], tag_to_content['AF_TGP'],
                                        tag_to_content['CLNDN'], tag_to_content['CLNHGVS'], tag_to_content['CLNSIG'],
                                        so[3:], molecular_consequence]
                    idx = 0
                    for ele in join_list:
                        if ele == '':
                            join_list[idx] = 'N'
                        idx += 1
                    new_row = ','.join(join_list)+'\n'
                    clinvar_table.write(new_row)
        else:
            clinvar_non_genename.write(row)
clinvar_table.close()
clinvar_non_genename.close()

gene_to_info = {}
for row in open('/home/administrator/database/mutation_database/Homo_sapiens.GRCh37.75.gtf-gene_info.gtf','r'):
    chr, start, end, strand, ensg, genename = row.strip().split('\t')
    gene_to_info[genename] = ','.join([chr, start, end, ensg]) + ','

# final result
output = open('/home/administrator/database/mutation_database/clinvar/clinvar_breast_93genes_20180729.csv', 'w')
output.write('genename,chr,start,end,geneid,pos,ref,alt,alleleid,rs,af_esp,af_exac,af_tgp,clngn,clnhgvs,clnsig\n')

#match_gene = []
target_genes = []
for row in open('/home/administrator/database/mutation_database/qiagen_breast_cancer_panel_gene_list.txt', 'r'):
    target_genes.append(row.strip())

# exist_genes = set()
for row in open('/home/administrator/database/mutation_database/clinvar/clinvar_table.csv', 'r'):
    genename = row.split(',')[0]
    info = ','.join(row.strip().split(',')[1:])
    if genename in target_genes:
        infos = genename + ',' + gene_to_info[genename] + info + '\n'
        output.write(infos)
#         exist_genes.add(genename)
output.close()