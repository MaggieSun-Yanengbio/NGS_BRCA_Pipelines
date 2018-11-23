import pymysql

conn = pymysql.connect(
    host='localhost',
    user='root',
    passwd='admin',
    db='annotation'
)
cur = conn.cursor()
cur.execute('create table if not exists mutations(GENE varchar(100), CHR varchar(100), POSITION varchar(100), '
            'REF varchar(1000), ALT varchar(1000), HGVS varchar(1000) primary key, '
            'RS_DBSNP varchar(100), RS_1000GENOMES varchar(100), CLNDN_CLINVAR varchar(1000), '
            'AF_1000GENOMES varchar(100), EAS_AF_1000GENOMES varchar(100), AMR_AF_1000GENOMES varchar(100), AFR_AF_1000GENOMES varchar(100), EUR_AF_1000GENOMES varchar(100), SAS_AF_1000GENOMES varchar(100), AF_EXAC varchar(100), AN_EAS_EXAC varchar(100), AF_CNGB varchar(100), AF_LOCAL varchar(100), '
            'COSMID_COSMIC varchar(100), FEATURE_ID_COSMIC varchar(100), CDS_LENGTH_COSMIC varchar(100), CDS_LENGTH_EVS varchar(1000), HGVS_C_COSMIC varchar(100), HGVS_C_EVS varchar(1000), HGVS_P_COSMIC varchar(100), HGVS_P_EVS varchar(1000), '
            'MUTATION_DESCRIPTION_CLINVAR varchar(100),MUTATION_DESCRIPTION_COSMIC varchar(100), MUTATION_DESCRIPTION_EVS varchar(1000), '
            'SIGNIFICANCE_CLINVAR varchar(100), SIGNIFICANCE_UTAHDB varchar(100), FATHMM_PREDICT_COSMIC varchar(100), FATHMM_SCORE_COSMIC varchar(100), POLYPHEN2_PREDICT_EVS varchar(1000), '
            'PMID_COSMIC varchar(100), REFERENCE1_UTAHDB varchar(100), REFERENCE2_UTAHDB varchar(100), PM_DBSNP varchar(10), '
            'SO_CLINVAR varchar(100), GRCH38_POSITION_EVS varchar(100));')

cur.close()
conn.commit()
conn.close()
