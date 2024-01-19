def _params_get_ensembl(wildcards):
  
    if GENCODE_RELEASE.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'
    return organism

rule gft_exonBED:
    input:
        'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
    output:
        'resources/external/gencode_{release}/{genome}.exons.bed'
    log:
        'logs/external_data/gencode_{release}_{genome}_gft_exonBED.log'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_biotype
        """
            zgrep -P "\texon\t" {input} | cut -f1,4,5,7,9 | \
            sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$6,".",$4 }}' | \
            sort -k1,1 -k2,2n > {output} 
        """



rule gft_transcriptBED:
    input:
        'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
    output:
        'resources/external/gencode_{release}/{genome}.transcripts.bed'
    log:
        'logs/external_data/gencode_{release}_{genome}_gft_transcriptBED.log'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
            zgrep -P "\ttranscript\t" {input} | cut -f1,4,5,7,9 | \
            sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$6,".",$4 }}' | \
            sort -k1,1 -k2,2n > {output} 
        """


rule gft_geneBED:
    input:
        'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
    output:
        'resources/external/gencode_{release}/{genome}.geneset.tsv'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
            zgrep -P "\tgene\t" {input} | cut -f1,4,5,7,9 | \
            sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$6,$4,$3-$2+1,$10,$8 }}' | \
            sort -k1,1 -k2,2n > {output} 
        """