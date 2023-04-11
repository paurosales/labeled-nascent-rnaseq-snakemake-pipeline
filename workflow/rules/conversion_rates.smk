# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('external/gencode/Mus_musculus.GRC{genome}.primary_assembly.genome.fa', genome=GENOME)

def _input_bedFile(wildcards):
    return expand('external/gencode/Mus_musculus.GRC{genome}.transcripts.bed', genome=GENOME)



rule alley_rates:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        'results/{proj}/alleyoop/rates/{label}_mapped_filtered_overallrates.csv'
    params:
        outdir = 'results/{proj}/alleyoop/rates',
        extra_params = ''  # [-mq <MQ cutoff>]
    resources:
        mem_mb = 2000
    threads: 6 
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            alleyoop rates -o {params.outdir} -r {input.ref_genome} -t {threads} {params.extra_params} {input.filteredBAM}
        """

rule alley_utrrates:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        'results/{proj}/alleyoop/utrrates/{label}_mapped_filtered_mutationrates_utr.csv',
        'results/{proj}/alleyoop/utrrates/{label}_mapped_filtered_mutationrates_utr.pdf'
    params:
        outdir = 'results/{proj}/alleyoop/utrrates',
        BED_file = _input_bedFile,
        extra_params = '-l 100 -m'  # [-mq <MQ cutoff> default = 27] [-m]
    resources:
        mem_mb = 12000
    threads: 6 
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            alleyoop utrrates -o {params.outdir} -r {input.ref_genome} -b {params.BED_file} -t {threads} {params.extra_params} {input.filteredBAM}
        """