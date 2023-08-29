# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME)

def _input_bedFile(wildcards):
    sample_type = wildcards['sample_type']
    treatment = wildcards['treatment']
    bio_rep = wildcards['bio_rep']

    seq_mode = SAMPLES.loc[( 
                          sample_type, 
                          treatment,
                          bio_rep), 
                        'Seq_mode']

    if seq_mode == 'mRNA':
        bed_file = expand('resources/external/gencode_{release}/{genome}.transcripts.bed', release=GENCODE_RELEASE, genome=GENOME)
    
    elif seq_mode == '3prime':
        bed_file = expand('resources/external/gencode_{release}/{genome}.3UTR.bed', release=GENCODE_RELEASE, genome=GENOME)

    else:
        bed_file = expand('resources/external/gencode_{release}/{genome}.transcripts.bed', release=GENCODE_RELEASE, genome=GENOME)

    return bed_file


def _params_snpe_eval(wildcards):

    sample_type = wildcards['sample_type']
    treatment = wildcards['treatment']
    bio_rep = wildcards['bio_rep']

    seq_len = SAMPLES.loc[( 
                          sample_type, 
                          treatment,
                          bio_rep), 
                        'Seq_length']

    seq_len += 10

    return seq_len


rule alley_snpeval:
    input:
        snpVCF = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_snp.vcf',
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome,
        BED_file = _input_bedFile
    output:
        'results/quality_control/snp_eval/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_SNPeval.csv',
        'results/quality_control/snp_eval/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_SNPeval.pdf'
    params:
        outdir = 'results/quality_control/snp_eval/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        snp_dir = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        var_cov = config['SLAM']['VAR_MIN_COVERAGE'],
        var_fract = config['SLAM']['VAR_FRACTION'],
        max_len = _params_snpe_eval,
        min_qual = config['SLAM']['MIN_QUAL']
    resources:
        mem_mb = 2000
    threads: 6 
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop snpeval -o {params.outdir} -s {params.snp_dir} \
            -r {input.ref_genome} -b {input.BED_file} \
            -c {params.var_cov} -f {params.var_fract} \
            -l {params.max_len} -q {params.min_qual} -m \
            -t {threads} {input.filteredBAM}
        """

rule alley_tcperreadpos:
    input: 
        snpVCF = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_snp.vcf',
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        'results/quality_control/tc_position/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_tcperreadpos.csv',
        'results/quality_control/tc_position/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_tcperreadpos.pdf'
    params:
        outdir = 'results/quality_control/tc_position/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        snp_dir = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        min_qual = config['SLAM']['MIN_QUAL'],
        max_len = _params_snpe_eval,
    resources:
        mem_mb = 12000
    threads: 6 
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop tcperreadpos -o {params.outdir} \
            -r {input.ref_genome} -s {params.snp_dir} \
            -l {params.max_len} -mq {params.min_qual} \
            -t {threads} {input.filteredBAM}
        """


rule alley_tcperutrpos:
    input: 
        snpVCF = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_snp.vcf',
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome,
        BED_file = _input_bedFile
    output:
        'results/quality_control/tc_position/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_tcperutr.csv',
        'results/quality_control/tc_position/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_tcperutr.pdf'
    params:
        outdir = 'results/quality_control/tc_position/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        snp_dir = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        min_qual = config['SLAM']['MIN_QUAL'],
        max_len = _params_snpe_eval,
    resources:
        mem_mb = 12000
    threads: 6 
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop tcperutrpos -o {params.outdir} \
            -r {input.ref_genome} -b {input.BED_file} -s {params.snp_dir} \
            -l {params.max_len} -mq {params.min_qual} \
            -t {threads} {input.filteredBAM}
        """