# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{realease}/{genome}.genome.fa', realease=GENCODE_REALEASE, genome=GENOME)

def _input_bedFile(wildcards):
    return expand('resources/external//gencode_{realease}/{genome}.transcripts.bed', genome=GENOME)

def _params_for_ngm_mapPE1(wildcards):
    label = wildcards.label
    row = SAMPLE_MANIFEST.loc[label]
    sampleID = row['Sample_ID']
    return sampleID


def _params_for_ngm_mapPE2(wildcards):
    label = wildcards.label
    row = SAMPLE_MANIFEST.loc[label]
    sampleName = row['Sample_Name']
    return sampleName


rule ngm_mapPE:
    input: 
        ref_genome = _input_refGenome,
        fq_trimmed_1 = 'results/trimmed_fastq/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq_trimmed_2 = 'results/trimmed_fastq/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        mapBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}.bam'
    params:
        label = '{sample_type}_{treatment}_Bio-rep_{bio_rep}'
    resources:
        mem_mb = 20000 # 8 GB for human (4 GB for genomes < 500 MBp) according to https://github.com/Cibiv/NextGenMap/wiki
    threads: 36
    log:
        'logs/{sample_type}_{treatment}_Bio-rep_{bio_rep}_ngm.log'
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            ngm -b -r {input.ref_genome}\
            -1 {input.fq_trimmed_1}\
            -2 {input.fq_trimmed_2}\
            -t {threads} --no-progress --slam-seq 2 -5 12\
            --rg-id {params.label} --rg-sm {params.label}:NA:-1\
            -o {output.mapBAM} 2>&1 > {log}
        """


rule slam_filter:
    input: 
        BED_file = _input_bedFile,
        mapBAM = 'results/slamdunk/map/{sample_type}_{treatment}_Bio-rep_{bio_rep}.bam'
    output:
        outdir = directory('results/bam_files'),
        filteredBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam',
        indexBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam.bai'
    params:
        min_qual = config['SLAM']['MIN_MAP_QUALITY'],
        min_ident = config['SLAM']['MIN_IDENTITY'],
        max_mismatch = config['SLAM']['MAX_MISMATCH']
    resources:
        mem_mb = 8000
    threads: 30
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            slamdunk filter -o {output.outdir} -b {input.BED_file}\
            -mq {params.min_qual} -mi {params.min_ident} -nm {params.max_mismatch}\
            -t {threads} {input.mapBAM}
        """