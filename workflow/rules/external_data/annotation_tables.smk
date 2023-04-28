def _params_get_ensembl(wildcards):
  
    if GENCODE_RELEASE.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'
    return organism


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
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$8,".",$4,$6 }}' | \
            sort -k1,1 -k2,2n > {output} 
        """


# zgrep -P "\ttranscript\t" GRCm38.annotation.gtf.gz | cut -f1,4,5,7,9 | \
# sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
# awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2,$3,$8,".",$4,$6 }' | \
# sort -k1,1 -k2,2n > GRCm38.transcripts.bed

# chr2    HAVANA  transcript      87517032        87517500        .       -       .       
# gene_id "ENSMUSG00000081762.1"; 
# transcript_id "ENSMUST00000122437.1"; # 8
# gene_type "processed_pseudogene"; 
# gene_name "Gm13738"; 
# transcript_type "processed_pseudogene"; # 14
# transcript_name "Gm13738-201"; 
# level 1; 
# transcript_support_level "NA"; 
# mgi_id "MGI:3649860"; 
# ont "PGO:0000004"; 
# tag "pseudo_consens"; 
# tag "basic"; 
# havana_gene "OTTMUSG00000013807.1"; 
# havana_transcript "OTTMUST00000033014.1";

rule gft_geneBED:
    input:
        'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
    output:
        'resources/external/gencode_{release}/{genome}.genes.bed'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
            zgrep -P "\tgene\t" {input} | cut -f1,4,5,7,9 | \
            sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$6,".",$4,$3-$2+1,$10,$8 }}' | \
            sort -k1,1 -k2,2n > {output} 
        """


# chr16	HAVANA	gene	8513429	8621568	.	+	.	
# gene_id "ENSMUSG00000057880.12"; 
# gene_type "protein_coding"; 
# gene_name "Abat"; 
# level 2; 
# mgi_id "MGI:2443582"; 
# havana_gene "OTTMUSG00000025937.6";


# zgrep -P "\tgene\t" GRCm38.annotation.gtf.gz | cut -f1,4,5,7,9 | \
# sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
# awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2,$3,$6,".",$4,$3-$2+1,$10,$8 }' | \
# sort -k1,1 -k2,2n > AUX2.gene.bed
