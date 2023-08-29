# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa.gz', release=GENCODE_RELEASE, genome=GENOME)

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
        bed_file = expand('resources/external/gencode_{release}/{genome}.exons.bed', release=GENCODE_RELEASE, genome=GENOME)
    
    elif seq_mode == '3prime':
        bed_file = expand('resources/external/gencode_{release}/{genome}.3UTR.bed', release=GENCODE_RELEASE, genome=GENOME)

    else:
        bed_file = expand('resources/external/gencode_{release}/{genome}.transcripts.bed', release=GENCODE_RELEASE, genome=GENOME)

    return bed_file


rule ngm_mapPE:
    input: 
        ref_genome = _input_refGenome,
        fq_trimmed_1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq_trimmed_2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        mapBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.bam'
    params:
        label = '{sample_type}_{treatment}_Bio-rep_{bio_rep}'
    resources:
        mem_mb = 20000 # 8 GB for human (4 GB for genomes < 500 MBp) according to https://github.com/Cibiv/NextGenMap/wiki
    threads: 30
    log:
        'logs/mapping/{sample_type}_{treatment}_Bio-rep_{bio_rep}_ngm.log'
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            ngm -b -r {input.ref_genome}\
            -1 {input.fq_trimmed_1}\
            -2 {input.fq_trimmed_2}\
            -t {threads} --no-progress --slam-seq 2 \
            --rg-id {params.label} --rg-sm {params.label}:NA:-1 \
            -o {output.mapBAM} > {log}  2>&1
        """


rule slam_filter:
    input: 
        BED_file = _input_bedFile,
        mapBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.bam'
    output:
        filteredBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam',
        indexBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam.bai'
    params:
        outdir = directory('results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        min_qual = config['SLAM']['MIN_MAP_QUALITY'],
        min_ident = config['SLAM']['MIN_IDENTITY'],
        max_mismatch = config['SLAM']['MAX_MISMATCH']
    resources:
        mem_mb = 12000
    threads: 20
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            slamdunk filter -o {params.outdir} -b {input.BED_file}\
            -mq {params.min_qual} -mi {params.min_ident} -nm {params.max_mismatch}\
            -t {threads} {input.mapBAM}
        """