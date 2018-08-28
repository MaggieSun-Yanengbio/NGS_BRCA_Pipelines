#cosmic数据整理
filein = open('/home/administrator/database/mutation_database/COSMIC/CosmicMutantExport_20180820.tsv','r')
output = open('/home/administrator/database/mutation_database/COSMIC/cosmic_breast_93genes_20180820.csv','w')

target_genes = []
for gene in open('/home/administrator/database/mutation_database/qiagen_breast_cancer_panel_gene_list.txt','r'):
    target_genes.append(gene.strip())

filein.readline()
for row in filein.readlines():
    genename, accessing_number, cds_length, hgnc_id, sample_name, sample_id, tumor_id, primary_site, site_sub1, \
    site_sub2 , site_sub3, primary_histology, hist_sub1, hist_sub2, hist_sub3, genome_wide_screen, mutation_id, \
    mutation_cds, mutation_aa, mutation_description, mutation_zygosity, log, grch, pos, mutation_strand, snp, \
    resistance, prediction, score, mutation_somatic_status, pubmed, id_study, sample_type, tumor_origin, age \
    = row.split('\t')
    if 'ENST' in genename:
        genename = genename[:genename.find('_')]
    if genename == 'MRE11A':
        genename = 'MRE11'
    if genename in target_genes:
        if hgnc_id.strip() == '':
            hgnc_id = 'NS'
        if '?' in mutation_cds:
            mutation_cds = 'NS'
        if '?' in mutation_aa:
            mutation_aa = 'NS'
        if 'Unknown' in mutation_description:
            mutation_description = 'NS'
        if mutation_zygosity.strip() == '':
            mutation_zygosity = 'NS'
        if grch.strip() == '':
            grch = 'NS'
        if pos.strip() == '':
            chr = start = end = 'NS'
        else:
            chr = pos[:pos.find(':')]
            start = pos[pos.find(':')+1:pos.find('-')]
            end = pos[pos.find('-')+1:]
        if '-' in resistance:
            resistance = 'NS'
        if mutation_strand.strip() == '':
            mutation_strand = 'NS'
        if snp.strip() == '':
            snp = 'NS'
        if prediction.strip() == '':
            prediction = 'NS'
        if score.strip() == '':
            score = '-1'
        if pubmed.strip() == '':
            pubmed = 'NS'
        if id_study.strip() == '':
            id_study = 'NS'
        if age.strip() == '':
            age = '-1\n'
        new_row = ','.join([genename, accessing_number, cds_length, hgnc_id, sample_name, sample_id, tumor_id,
                            primary_site, site_sub1, site_sub2 , site_sub3, primary_histology, hist_sub1,
                            hist_sub2, hist_sub3, genome_wide_screen, mutation_id, mutation_cds, mutation_aa,
                            mutation_description, mutation_zygosity, log, grch, chr, start, end, mutation_strand,
                            snp, resistance, prediction, score, mutation_somatic_status, pubmed, id_study, sample_type,
                            tumor_origin, age])
output.write(new_row)
