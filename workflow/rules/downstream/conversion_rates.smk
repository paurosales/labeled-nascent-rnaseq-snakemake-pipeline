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


rule alley_rates:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_overallrates.csv',
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_overallrates.pdf'
    params:
        outdir = 'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        min_qual = config['SLAM']['MIN_QUAL']
    resources:
        mem_mb = 2000
    threads: 6 
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop rates -o {params.outdir} -r {input.ref_genome} \
            -t {threads} -mq {params.min_qual} \
            {input.filteredBAM}
        """

rule alley_utrrates:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_mutationrates_utr.csv',
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_mutationrates_utr.pdf'
    params:
        outdir = 'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        BED_file = _input_bedFile,
        min_qual = config['SLAM']['MIN_QUAL'],
        max_len = _params_snpe_eval,
    resources:
        mem_mb = 12000
    threads: 6 
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop utrrates -o {params.outdir} -r {input.ref_genome} \
            -b {params.BED_file} -t {threads} \
            -l {params.max_len} -m -mq {params.min_qual} \
            {input.filteredBAM}
        """