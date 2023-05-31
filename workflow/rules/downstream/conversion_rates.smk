# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME)

def _input_bedFile(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.transcripts.bed', release=GENCODE_RELEASE, genome=GENOME)



rule alley_rates:
    input: 
        filteredBAM = rules.slam_filter.output.filteredBAM,
        ref_genome = _input_refGenome
    output:
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_overallrates.csv',
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_overallrates.pdf'
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
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_mutationrates_utr.csv',
        'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_mutationrates_utr.pdf'
    params:
        outdir = 'results/conversion_rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        BED_file = _input_bedFile,
        min_qual = config['SLAM']['MIN_QUAL'],
        max_len = config['SLAM']['MAX_LENGTH']
    resources:
        mem_mb = 12000
    threads: 6 
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop utrrates -o {params.outdir} -r {input.ref_genome} \
            -b {params.BED_file} -t {threads} \
            -l {params.max_len} -m {params.min_qual} \
            {input.filteredBAM}
        """