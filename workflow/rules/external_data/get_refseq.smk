def _params_get_seq(wildcards):
  
    if wildcards.realease.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'
    return organism
    # return dict(organism=organism, link_handle=link_handle) # use params[0][organism] o shell



rule get_genome:
    output:
        'resources/external/gencode_{realease}/GRC{genome}.genome.fa.gz'
    params:
        _params_get_seq
    log:
        'logs/external_data/get_gencode_{realease}_{genome}_genome.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params.organism}/release_{wildcards.realease}/{wildcards.genome}.primary_assembly.genome.fa.gz -O  {output}
        """


rule unzip_genome:
    input:
        rules.get_genome.output
    output:
        temp('resources/external/gencode_{realease}/GRC{genome}.genome.fa')
    shell:
        """
             gzip -dc {input} > {output}
        """

rule get_annotation:
    output:
        'resources/external/gencode_{realease}/GRC{genome}.annotation.gtf.gz'
    params:
        organism = _params_get_seq
    log:
        'logs/external_data/get_gencode_{realease}_{genome}_annotation.log'
    threads: 1
    resources:
        http = 1
    shell:
        """
            wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params.organism}/release_{wildcards.realease}/gencode.v{wildcards.realease}.annotation.gft.gz -O  {output}
        """


rule unzip_gft:
    input:
        rules.get_genome.output
    output:
        temp('resources/external/gencode_{realease}/GRC{genome}.annotation.gft')
    shell:
        """
             gzip -dc {input} > {output}
        """