rule alley_separator:
    input:
        snpVCF = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_snp.vcf',
        filteredBAM = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam',
        ref_genome = _input_refGenome
    output:
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_backgroundReads.bam',
        'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_TCReads.bam'
    params:
        outdir = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        snp_dir = 'results/variant_call/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        conv_th = config['SLAM']['CONVERSION_THRESHOLD'],
        min_qual = config['SLAM']['MIN_QUAL']
    resources:
        mem_mb = 12000
    threads: 32
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop read-separator -o {params.outdir} \
            -s {params.snp_dir} -r {input.ref_genome} \
            -c {params.conv_th} -q {params.min_qual}\
            -t {threads} {input.filteredBAM}
        """
