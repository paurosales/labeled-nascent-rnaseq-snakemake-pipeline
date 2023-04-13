def _deseq2_input(wildcards):
    return expand('results/{proj}/slamdunk/count/{label}_mapped_filtered_tcount.tsv', proj=PROJECT, label=LABELS)


rule deseq2:
    input:
        countFolder = _deseq2_input,
        sample_man = 'config/sample_manifest.tsv'
    output:
        fig_outdir = directory('results/figures/dif_expr'), # CHANGE FOR PATH INCLUDING PLOT.PDF
        data_outdir = directory('results/dif_expr')
        # PCA, MA plots and DESeq2.txt (contrast? folder depends on the design)
    params:
        subgroup_filter = config['deseq2']['subgroup_filter'],
        pval_th = config['deseq2']['pval_threshold']
    resources:
        mem_mb = 4000
    log: 
        'logs/{proj}/downstream/deseq2.log'
    conda:
        '../../envs/deseq2.yaml'
    script:
        "../../scripts/downstream/DE_analysis.R"
        