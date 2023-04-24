def _params_get_ensembl(wildcards):
  
    if GENCODE_RELEASE.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'
    return organism

rule gtf2bed:
    input:
        'resources/external/gencode_{release}/{genome}.annotation.gtf'
    output:
        'resources/external/gencode_{release}/{genome}.transcripts.bed'
    log:
        'logs/external_data/gencode_{release}_{genome}_gtf2bed.log'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
        zgrep -P "\ttranscript\t" {input} | cut -f1,4,5,7,9 | \
        sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
        awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2-1,$3,$10,".",$4,$14,$12,$18 }}' | \
        sort -k1,1 -k2,2n > {output} > {log}  2>&1
        """


# zgrep -P "\ttranscript\t" Mus_musculus.m39.108.gtf.gz | cut -f1,4,5,7,9 | \
# sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
# awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,($10)"."($12),".",$4,$14,$18 }' | \
# sort -k1,1 -k2,2n > Mus_musculus.m39_108.transcripts.bed

# TRANSCRIPT ATTRIBUTES
# gene_id "ENSMUSG00000095742"; 
# gene_version "2"; 
# transcript_id "ENSMUST00000179436"; 
# transcript_version "2"; 
# exon_number "7"; 
# gene_source "ensembl"; 
# gene_biotype "protein_coding"; 
# transcript_source "ensembl"; 
# transcript_biotype "protein_coding"; #18
# exon_id "ENSMUSE00001007635"; 
# exon_version "2"; 
# tag "basic"; 
# tag "Ensembl_canonical"; 
# transcript_support_level "5 (assigned to previous version 1)"
 
# CHANGEEEEEE
rule get_ensembl_geneset:
    output:
        genesetTSV = 'resources/external/ensembl/{genome}.geneset.tsv',
        transcriptsetTSV = 'resources/external/ensembl/{genome}.transcriptset.tsv'
    log:
        'logs/external_data/{genome}_ensembl_geneset.log'
    params: 
        ensembl_version = config['ENSEMBL']['VERSION'],
        organism = _params_get_ensembl
    threads: 6
    conda:
        '../../envs/downstream/biomart.yaml'
    resources:
        mem_mb = 4000
    script:
        "../../scripts/external_data/ensembl_geneset.R"