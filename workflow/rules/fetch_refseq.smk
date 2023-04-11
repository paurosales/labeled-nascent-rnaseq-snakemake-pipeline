rule get_refgenome:
    output:
        'external/gencode/Mus_musculus.GRC{genome}.primary_assembly.genome.fa'
    params:
        genome_build =  config['genome_build'],
        gencode_version =  config['gencode_version']
    shell:
        """
        wget -qO http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{gencode_version}/GRC{genome_build}.primary_assembly.genome.fa.gz | gunzip -c > {output} 
        """


rule get_annotation:
    output:
        'external/gencode/Mus_musculus.GRC{genome}.annotation.gtf'
    params:
        genome_build =  config['genome_build'],
        gencode_version =  config['gencode_version']
    shell:
        """
        wget -qO http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{gencode_version}/gencode.vM{gencode_version}.annotation.gtf.gz | gunzip -c > {output} 
        """


rule gtf2bed:
    input:
        'external/gencode/Mus_musculus.GRC{genome}.annotation.gtf'
    output:
        'external/gencode/Mus_musculus.GRC{genome}.transcripts.bed'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
        zgrep -P "\ttranscript\t" {input} | cut -f1,4,5,7,9 | \
        sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
        awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$10,".",$4,$14,$12,$18 }' | \
        sort -k1,1 -k2,2n > {output}
        """


# zgrep -P "\ttranscript\t" Mus_musculus.GRCm39.108.gtf.gz | cut -f1,4,5,7,9 | \
# sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
# awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,($10)"."($12),".",$4,$14,$18 }' | \
# sort -k1,1 -k2,2n > Mus_musculus.GRCm39_108.transcripts.bed

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
 