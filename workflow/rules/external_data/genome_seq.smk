rule get_refgenome:
    output:
        'resources/external/gencode/Mus_musculus.GRC{genome}.primary_assembly.genome.fa'
    params:
        genome_build =  config['genome_build'],
        gencode_version =  config['gencode_version']
    shell:
        """
        wget -qO http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{gencode_version}/GRC{genome_build}.primary_assembly.genome.fa.gz | gunzip -c > {output} 
        """


rule get_annotation:
    output:
        'resources/external/gencode/Mus_musculus.GRC{genome}.annotation.gtf'
    params:
        genome_build =  config['genome_build'],
        gencode_version =  config['gencode_version']
    shell:
        """
        wget -qO http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{gencode_version}/gencode.vM{gencode_version}.annotation.gtf.gz | gunzip -c > {output} 
        """