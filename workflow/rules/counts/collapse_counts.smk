def _input_genset(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)


rule alley_collapse:
    input: 
        'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_tcount.tsv'
    output:
        counts = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filtered_tcount_collapsed.csv'
    params:
        outdir = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}'
    resources:
        mem_mb = 12000
    threads: 16
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop collapse -o {params.outdir} -t {threads} {input}
        """


rule countmtx_TotalReads:
    input: 
        raw_cts = TARGETS['slam_count'],
        collapsed_cts = TARGETS['alley_collapse'],
        geneset_TSV = _input_genset
    output:
        'results/counts/count_matrix/TotalReadCount_mtx.tsv'
    params:
        eval_expression = 'ReadCount'
    resources:
        mem_mb = 4000
    threads: 4
    log:
        'logs/counts/count_matrix/TotalReadCount_merge.log'
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    script:
        "../../scripts/counts/collapse_counts.R"



rule countmtx_TcReads:
    input: 
        raw_cts = TARGETS['slam_count'],
        collapsed_cts = TARGETS['alley_collapse'],
        geneset_TSV = _input_genset
    output:
        'results/counts/count_matrix/TcReadCount_mtx.tsv'
    params:
        eval_expression = 'TcReadCount'
    resources:
        mem_mb = 4000
    threads: 4
    log:
        'logs/counts/count_matrix/TcReadCount_merge.log'
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    script:
        "../../scripts/counts/collapse_counts.R"


rule countmtx_Tc_and_Total:
    input: 
        raw_cts = TARGETS['slam_count'],
        collapsed_cts = TARGETS['alley_collapse'],
        geneset_TSV = _input_genset
    output:
        'results/counts/count_matrix/TcReadCount_to_TotalReadCount_mtx.tsv'
    params:
        eval_expression = 'TcReadCount/ReadCount'
    resources:
        mem_mb = 4000
    threads: 4
    log:
        'logs/counts/count_matrix/TcReadCount_to_TotalReadCount_merge.log'
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    script:
        "../../scripts/counts/collapse_counts.R"

