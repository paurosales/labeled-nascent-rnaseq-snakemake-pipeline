# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('external/gencode/Mus_musculus.GRC{genome}.primary_assembly.genome.fa', genome=GENOME)

def _input_bedFile(wildcards):
    return expand('external/gencode/Mus_musculus.GRC{genome}.transcripts.bed', genome=GENOME)


rule slam_snp:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        snpVCF = 'results/{proj}/slamdunk/snp/{label}_mapped_filtered_snp.vcf'
    params:
        outdir = 'results/{proj}/slamdunk/snp',
        variant_fraction = config['slam']['variant_fraction'],
        extra_params = ''  # [-c <coverage cutoff>]
    resources:
        mem_mb = 6000
    threads: 18
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            slamdunk snp  -o {params.outdir} -r {input.ref_genome} -t {threads} -f {params.variant_fraction} {params.extra_params} {input.filteredBAM}
        """


rule slam_count:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        snpVCF = rules.slam_snp.output.snpVCF,
        ref_genome = _input_refGenome
    output:
        'results/{proj}/slamdunk/count/{label}_mapped_filtered_tcount.tsv'
    params:
        outdir = 'results/{proj}/slamdunk/count',
        snp_dir = 'results/{proj}/slamdunk/snp',
        BED_file = _input_bedFile,
        max_len = config['slam']['max_len'],
        conv_threshold = config['slam']['conversion_threshold'],
        min_qual = config['slam']['min_qual']
    resources:
        mem_mb = 8000
    threads: 30
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            slamdunk count -o {params.outdir} -s {params.snp_dir} -r {input.ref_genome} -b {params.BED_file} -t {threads}\
            -l {params.max_len} -c {params.conv_threshold} -q {params.min_qual} {input.filteredBAM}
        """
