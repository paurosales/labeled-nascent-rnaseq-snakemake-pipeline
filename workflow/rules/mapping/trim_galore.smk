

rule trim_galore:
    input:
        fq_1 = 'resources/raw_data/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz',
        fq_2 = 'resources/raw_data/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz'
    output:
        out_dir = directory('results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        fastqc_dir = directory('results/quality_control/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        report1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz_trimming_report.txt',
        report2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz_trimming_report.txt',
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    threads: 16
    resources:
        mem_mb = 5000
    conda:
        '../../envs/raw_processing/trim-galore.yaml'
    log: 
        'logs/fastq/trim_galore/{sample_type}_{treatment}_Bio-rep_{bio_rep}.log'
        # str(LOG_DIR / 'trim_galore' / f'{sample_type}_{treatment}_Bio-rep_{bio_rep}.log')
    params:
        basename = '{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        clip_5 = config['TRIM_GALORE']['CLIP_5'],
        clip_3 = config['TRIM_GALORE']['CLIP_3'],
        extra = '--2colour 20'  # SOLVE PARAMS PARSING
    shell:
        """
            mkdir -p {output.fastqc_dir} && \
            trim_galore {params.extra} -j {threads} \
            --basename {params.basename} \
            --clip_R1 {params.clip_5} --clip_R2 {params.clip_5} \
            --three_prime_clip_R1 {params.clip_3} --three_prime_clip_R1 {params.clip_3} \
            --gzip --paired --fastqc_args "--outdir {output.fastqc_dir}"\
            -o {output.out_dir} {input} > {log}  2>&1
        """

