# Change of files handle
def _input_for_trim_fastq(wildcards):
    label = wildcards.label
    fq1, fq2 = raw_fastq_identifiers_to_label(label)
    return {'fastq1': fq1, 'fastq2': fq2}

rule trimgalore:
    input:
        unpack(_input_for_trim_fastq)       # uses identifiers as handles and then changes to labels
    output:
        fq_trimmed_1 = 'results/trimmed_fastq/{label}_val_1.fq.gz',
        fq_trimmed_2 = 'results/trimmed_fastq/{label}_val_2.fq.gz'
    threads: 15
    resources:
        mem_mb = 8000
    conda:
        '../envs/trim.yaml'
    log: 
        'logs/trimmed_fastq/{label}_trim.log'
    params:
        outdir = 'results/trimmed_fastq',
        clip_5 = config['trim_galore']['clip_5'],
        clip_3 = config['trim_galore']['clip_3'],
        length = config['trim_galore']['length'],
        quality = config['trimgalore']['quality'], # [--2colour/--quality] depending on sequencer
        extra_params = config['trim_galore']['extra_params']
    shell:
        """
            trim_galore -j {threads} \
            --basename {wildcards.label} \
            --clip_R1 {params.clip_5} --clip_R2 {params.clip_5} \
            --three_prime_clip_R1 {params.clip_3} --three_prime_clip_R1 {params.clip_3} \
            --length {params.length} --fastqc --paired \
            {params.quality} \
            {params.extra_params} \
            -o {params.outdir} {input} > {log} 2>&1
        """
