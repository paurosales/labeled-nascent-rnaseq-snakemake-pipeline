def _deseq2_input(wildcards):
    return TARGETS['gene_counts']

rule deseq2:
    input:
        countFiles = _deseq2_input,
        sample_manifest = 'config/sample_manifest.tsv'
    output:
        fig_outdir = directory('results/figures/dif_expr'), # CHANGE FOR PATH INCLUDING PLOT.PDF
        data_outdir = directory('results/dif_expr'),
        ddsRDS = 'results/dif_expr/dds.Rds'
        # PCA, MA plots and DESeq2.txt (contrast? folder depends on the design)
    params:
        subgroup_filter = config['DESEQ2']['SUBSET_FILTER'],
        pval_th = config['DESEQ2']['PVAL_THRESHOLD']
    resources:
        mem_mb = 4000
    log: 
        'logs/downstream/deseq2.log'
    conda:
        '../../envs/downstream/deseq2.yaml'
    script:
        "../../scripts/downstream/DE_analysis.R"
        