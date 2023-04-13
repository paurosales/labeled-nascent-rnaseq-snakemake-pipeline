# Handle wildcards errors
def _input_for_TCcounts(wildcards):
    return expand('results/slamdunk/count/{label}_tcount_by_gene.tsv', proj=PROJECT, label=SAMPLE_MANIFEST.index)



# All-samples plots
rule all_TCcounts_plot:
    input:
        geneCounts = _input_for_TCcounts
    output:
        tcPercentPDF = 'results/figures/allSamples_TCReadsPercent_barplot.pdf'
    params:
        indir = 'results/slamdunk/count',
        genome_build = config['genome_build'],
        sample_manifest = config['samples']['manifest']
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/all_TCcounts_plot.log'
    conda:
        '../../envs/biomart.yaml'
    script:
        "../../scripts/TCreads_all_plots.R"
        

# Plots per sample
rule persample_TCcounts_plot:
    input:
        geneCountsRData = rules.counts_collapse_transcripts.output.geneCountsRData
    output:
        tcFreqPDF = 'results/figures/{label}/TCReadsFreq_barplot.pdf',
        tcFreqNoZeroPDF = 'results/figures/{label}/TCReadsFreqNoZero_barplot.pdf',
        exprCountsPDF = 'results/figures/{label}/TCReads_vs_Steady.pdf',
        exprCountsNoLimsPDF = 'results/figures/{label}/TCReads_vs_Steady_noLims.pdf',
        exprCpmPDF  = 'results/figures/{label}/TCReads_vs_Steady_CPM.pdf'
    params:
        genome_build = config['genome_build']
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/{label}_per-sample_plots.log'
    conda:
        '../../envs/biomart.yaml'
    script:
        "../../scripts/TCreads_per-sample_plots.R"

# All-samples plots (with gene length normalized counts)
rule all_normTCcounts_plot:
    input:
        geneCounts = _input_for_TCcounts
    output:
        tcPercentPDF = 'results/figures/norm_allSamples_TCReadsPercent_barplot.pdf'
    params:
        indir = 'results/slamdunk/count',
        sample_manifest = config['samples']['manifest']
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/norm_all_TCcounts_plot.log'
    conda:
        '../../envs/biomart.yaml'
    script:
        "../../scripts/figures/normTCreads_all_plots.R"
        

# Plots per sample (with gene length normalized counts)
rule persample_normTCcounts_plot:
    input:
        geneCountsRData = rules.counts_collapse_transcripts.output.geneCountsRData
    output:
        tcFreqPDF = 'results/figures/{label}/norm_TCReadsFreq_barplot.pdf',
        tcFreqNoZeroPDF = 'results/figures/{label}/norm_TCReadsFreqNoZero_barplot.pdf',
        exprCountsPDF = 'results/figures/{label}/norm_TCReads_vs_Steady.pdf',
        exprCountsNoLimsPDF = 'results/figures/{label}/norm_TCReads_vs_Steady_noLims.pdf',
        exprCpmPDF  = 'results/figures/{label}/norm_TCReads_vs_Steady_CPM.pdf'
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/{label}_norm_per-sample_plots.log'
    conda:
        '../../envs/biomart.yaml'
    script:
        "../../scripts/figures/normTCreads_per-sample_plots.R"
