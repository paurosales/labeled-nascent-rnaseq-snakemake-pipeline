rule bamCoverage_fwd_totalReads:
    input:
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam.bai'
    output:
        'results/bigwig_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_totalReads_coverage_fwd.bw'
    resources:
        mem_mb = 12000
    threads: 12
    conda:
        '../../envs/downstream/deeptools.yaml'
    log: 
        'logs/deeptools/{sample_type}_{treatment}_Bio-rep_{bio_rep}_totalReads_bamCoverage_fwd.log'
    params:
        bin_size = config['BAM_COVERAGE']['BIN_SIZE'],
        # smooth_length = config['BAM_COVERAGE']['SMOOTH_LENGTH'],
        mapq = config['BAM_COVERAGE']['MIN_MAPQ'],
        norm_reads = config['BAM_COVERAGE']['NORM_METHOD'],
        extra = config['BAM_COVERAGE']['EXTRA']
    shell:
        """
            bamCoverage --bam {input.bam} -o {output} \
            --filterRNAstrand forward \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --minMappingQuality {params.mapq} \
            --normalizeUsing {params.norm_reads} \
            {params.extra} > {log} 2>&1  
        """


rule bamCoverage_rev_totalReads:
    input:
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered.bam.bai'
    output:
        'results/bigwig_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_totalReads_coverage_rev.bw'
    resources:
        mem_mb = 12000
    threads: 12
    conda:
        '../../envs/downstream/deeptools.yaml'
    log: 
        'logs/deeptools/{sample_type}_{treatment}_Bio-rep_{bio_rep}_totalReads_bamCoverage_rev.log'
    params:
        bin_size = config['BAM_COVERAGE']['BIN_SIZE'],
        # smooth_length = config['BAM_COVERAGE']['SMOOTH_LENGTH'],
        mapq = config['BAM_COVERAGE']['MIN_MAPQ'],
        norm_reads = config['BAM_COVERAGE']['NORM_METHOD'],
        extra = config['BAM_COVERAGE']['EXTRA']
    shell:
        """
            bamCoverage --bam {input.bam} -o {output} \
            --filterRNAstrand reverse \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --minMappingQuality {params.mapq} \
            --normalizeUsing {params.norm_reads} \
            {params.extra} > {log} 2>&1  
        """

rule bamCoverage_fwd_backgroundReads:
    input:
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_backgroundReads.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_backgroundReads.bam.bai'
    output:
        'results/bigwig_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_backgroundReads_coverage_fwd.bw'
    resources:
        mem_mb = 12000
    threads: 12
    conda:
        '../../envs/downstream/deeptools.yaml'
    log: 
        'logs/deeptools/{sample_type}_{treatment}_Bio-rep_{bio_rep}_backgroundReads_bamCoverage_fwd.log'
    params:
        bin_size = config['BAM_COVERAGE']['BIN_SIZE'],
        # smooth_length = config['BAM_COVERAGE']['SMOOTH_LENGTH'],
        mapq = config['BAM_COVERAGE']['MIN_MAPQ'],
        norm_reads = config['BAM_COVERAGE']['NORM_METHOD'],
        extra = config['BAM_COVERAGE']['EXTRA']
    shell:
        """
            bamCoverage --bam {input.bam} -o {output} \
            --filterRNAstrand forward \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --minMappingQuality {params.mapq} \
            --normalizeUsing {params.norm_reads} \
            {params.extra} > {log} 2>&1  
        """


rule bamCoverage_rev_backgroundReads:
    input:
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_backgroundReads.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_backgroundReads.bam.bai'
    output:
        'results/bigwig_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_backgroundReads_coverage_rev.bw'
    resources:
        mem_mb = 12000
    threads: 12
    conda:
        '../../envs/downstream/deeptools.yaml'
    log: 
        'logs/deeptools/{sample_type}_{treatment}_Bio-rep_{bio_rep}_backgroundReads_bamCoverage_rev.log'
    params:
        bin_size = config['BAM_COVERAGE']['BIN_SIZE'],
        # smooth_length = config['BAM_COVERAGE']['SMOOTH_LENGTH'],
        mapq = config['BAM_COVERAGE']['MIN_MAPQ'],
        norm_reads = config['BAM_COVERAGE']['NORM_METHOD'],
        extra = config['BAM_COVERAGE']['EXTRA']
    shell:
        """
            bamCoverage --bam {input.bam} -o {output} \
            --filterRNAstrand reverse \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --minMappingQuality {params.mapq} \
            --normalizeUsing {params.norm_reads} \
            {params.extra} > {log} 2>&1  
        """

rule bamCoverage_fwd_TCReads:
    input:
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_TCReads.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_TCReads.bam.bai'
    output:
        'results/bigwig_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_TCReads_coverage_fwd.bw'
    resources:
        mem_mb = 12000
    threads: 12
    conda:
        '../../envs/downstream/deeptools.yaml'
    log: 
        'logs/deeptools/{sample_type}_{treatment}_Bio-rep_{bio_rep}_TCReads_bamCoverage_fwd.log'
    params:
        bin_size = config['BAM_COVERAGE']['BIN_SIZE'],
        # smooth_length = config['BAM_COVERAGE']['SMOOTH_LENGTH'],
        mapq = config['BAM_COVERAGE']['MIN_MAPQ'],
        norm_reads = config['BAM_COVERAGE']['NORM_METHOD'],
        extra = config['BAM_COVERAGE']['EXTRA']
    shell:
        """
            bamCoverage --bam {input.bam} -o {output} \
            --filterRNAstrand forward \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --minMappingQuality {params.mapq} \
            --normalizeUsing {params.norm_reads} \
            {params.extra} > {log} 2>&1  
        """


rule bamCoverage_rev_TCReads:
    input:
        bam = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_TCReads.bam',
        bai = 'results/bam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_TCReads.bam.bai'
    output:
        'results/bigwig_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_TCReads_coverage_rev.bw'
    resources:
        mem_mb = 12000
    threads: 12
    conda:
        '../../envs/downstream/deeptools.yaml'
    log: 
        'logs/deeptools/{sample_type}_{treatment}_Bio-rep_{bio_rep}_TCReads_bamCoverage_rev.log'
    params:
        bin_size = config['BAM_COVERAGE']['BIN_SIZE'],
        # smooth_length = config['BAM_COVERAGE']['SMOOTH_LENGTH'],
        mapq = config['BAM_COVERAGE']['MIN_MAPQ'],
        norm_reads = config['BAM_COVERAGE']['NORM_METHOD'],
        extra = config['BAM_COVERAGE']['EXTRA']
    shell:
        """
            bamCoverage --bam {input.bam} -o {output} \
            --filterRNAstrand reverse \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --minMappingQuality {params.mapq} \
            --normalizeUsing {params.norm_reads} \
            {params.extra} > {log} 2>&1  
        """
