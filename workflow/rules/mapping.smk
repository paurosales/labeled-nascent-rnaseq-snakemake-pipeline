# Handle wildcards errors
def _input_refGenome(wildcards):
    return expand('external/gencode/Mus_musculus.GRC{genome}.primary_assembly.genome.fa', genome=GENOME)

def _input_bedFile(wildcards):
    return expand('external/gencode/Mus_musculus.GRC{genome}.transcripts.bed', genome=GENOME)

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
        fq_trimmed_1 = 'results/{proj}/trimmed_fastq/{label}_val_1.fq.gz',
        fq_trimmed_2 = 'results/{proj}/trimmed_fastq/{label}_val_2.fq.gz'
    output:
        BAM = 'results/{proj}/slamdunk/map/{label}_mapped.bam'
        # SAM = 'results/{proj}/slamdunk/map/{label}_mapped.bam'
    params:
        sampleID = _params_for_ngm_mapPE1, 
        sampleName = _params_for_ngm_mapPE2
    resources:
        mem_mb = 20000 # 8 GB for human (4 GB for genomes < 500 MBp) according to https://github.com/Cibiv/NextGenMap/wiki
    threads: 36
    log:
        'logs/{proj}/{label}_ngm.log'
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            ngm -b -r {input.ref_genome}\
            -1 {input.fq_trimmed_1}\
            -2 {input.fq_trimmed_2}\
            -t {threads} --no-progress --slam-seq 2 -5 12\
            --rg-id {wildcards.label} --rg-sm {wildcards.label}:NA:-1\
            -o {output.BAM} 2>&1 > {log}
        """


rule slam_filter:
    input: 
        # rules.ngm_mapPE.output.BAM
        'results/{proj}/slamdunk/map/{label}_mapped.bam'
        # rules.remove_scaffolds.output.noScaffoldsBAM
    output:
        filteredBAM = 'results/{proj}/slamdunk/filter/{label}_mapped_filtered.bam',
        indexBAM = 'results/{proj}/slamdunk/filter/{label}_mapped_filtered.bam.bai'
    params:
        outdir = 'results/{proj}/slamdunk/filter',
        BED_file = _input_bedFile,
        extra_params = ''  # [-mq <MQ cutoff> (default: 2)] [-mi <identity cutoff> (default: 0.95)] [-nm <NM cutoff> (default: -1)]
    resources:
        mem_mb = 8000
    threads: 30
    conda:
        '../envs/slamdunk.yaml'
    shell:
        """
            slamdunk filter -o {params.outdir} -b {params.BED_file} -t {threads} {params.extra_params} {input}
        """
