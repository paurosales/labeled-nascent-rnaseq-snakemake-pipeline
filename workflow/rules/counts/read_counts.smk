# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{realease}/GRC{genome}.genome.fa', realease=GENCODE_REALEASE, genome=GENOME)

def _input_bedFile(wildcards):
    return expand('resources/external/gencode_{realease}/GRC{genome}.transcripts.bed', genome=GENOME)


rule slam_snp:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        snpVCF = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}.snp.vcf'
    params:
        outdir = 'results/variant_call',
        var_fract = config['SLAM']['VAR_FRACTION'],
        var_cov = config['SLAM']['VAR_MIN_COVERAGE']
    resources:
        mem_mb = 6000
    threads: 18
    conda:
        '../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            slamdunk snp  -o {params.outdir} -r {input.ref_genome} -t {threads} -f {params.car_fract} -c {params.var_cov} {input.filteredBAM}
        """


rule slam_count:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        snpVCF = rules.slam_snp.output.snpVCF,
        BED_file = _input_bedFile,
        ref_genome = _input_refGenome
    output:
        'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}.tcount.tsv'
    params:
        outdir = 'results/counts',
        snp_dir = 'results/variant_call',
        max_len = config['SLAM']['MAX_LENGTH'],
        conv_th = config['SLAM']['CONVERSION_THRESHOLD'],
        min_qual = config['SLAM']['MIN_QUAL']
    resources:
        mem_mb = 8000
    threads: 30
    conda:
        '../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            slamdunk count -o {params.outdir} -s {params.snp_dir} -r {input.ref_genome} -b {input.BED_file} -t {threads}\
            -l {params.max_len} -c {params.conv_th} -q {params.min_qual} {input.filteredBAM}
        """
