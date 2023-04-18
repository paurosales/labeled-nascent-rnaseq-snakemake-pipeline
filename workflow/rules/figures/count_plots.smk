# Handle wildcards errors CHANGEEEEE
def _input_for_TCcounts(wildcards):
    return expand('results/slamdunk/count/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.tsv', proj=PROJECT, label=SAMPLE_MANIFEST.index)

def _params_wildcard_genome(wildcards):
    return GENOME

# --------------------------- Using raw counts  --------------------------- #
# --- All samples plots
rule all_TCcounts_plot:
    input:
        geneCounts = _input_for_TCcounts
    output:
        tcPercentPDF = 'results/figures/allSamples_TCReadsPercent_barplot.pdf'
    params:
        # indir = 'results/slamdunk/count',
        genome_build = _params_wildcard_genome,
        sample_manifest = config['SAMPLE_MANIFEST']
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/all_TCcounts_plot.log'
    conda:
        '../../envs/downstream/biomart.yaml'
    script:
        "../../scripts/figures/TCreads_all_plots.R"
        

# --- Plots per sample
rule persample_TCcounts_plot:
    input:
        geneCountsRData = rules.counts_collapse_transcripts.output.geneCountsRData
    output:
        tcFreqPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/TCReadsFreq_barplot.pdf',
        tcFreqNoZeroPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/TCReadsFreqNoZero_barplot.pdf',
        exprCountsPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/TCReads_vs_Steady.pdf',
        exprCountsNoLimsPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/TCReads_vs_Steady_noLims.pdf',
        exprCpmPDF  = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/TCReads_vs_Steady_CPM.pdf'
    params:
        genome_build = _params_wildcard_genome
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}_per-sample_plots.log'
    conda:
        '../../envs/downstream/biomart.yaml'
    script:
        "../../scripts/figures/TCreads_per-sample_plots.R"

# --------------------------- Using normalized counts (gene size) --------------------------- #
# --- All samples plots
rule all_normTCcounts_plot:
    input:
        geneCounts = _input_for_TCcounts
    output:
        tcPercentPDF = 'results/figures/norm_allSamples_TCReadsPercent_barplot.pdf'
    params:
        indir = 'results/slamdunk/count',
        sample_manifest = config['SAMPLE_MANIFEST']
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/norm_all_TCcounts_plot.log'
    conda:
        '../../envs/downstream/biomart.yaml'
    script:
        "../../scripts/figures/normTCreads_all_plots.R"
        

# --- Plots per sample
rule persample_normTCcounts_plot:
    input:
        geneCountsRData = rules.counts_collapse_transcripts.output.geneCountsRData
    output:
        tcFreqPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/norm_TCReadsFreq_barplot.pdf',
        tcFreqNoZeroPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/norm_TCReadsFreqNoZero_barplot.pdf',
        exprCountsPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/norm_TCReads_vs_Steady.pdf',
        exprCountsNoLimsPDF = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/norm_TCReads_vs_Steady_noLims.pdf',
        exprCpmPDF  = 'results/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/norm_TCReads_vs_Steady_CPM.pdf'
    resources:
        mem_mb = 4000
    log: 
        'logs/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}_norm_per-sample_plots.log'
    conda:
        '../../envs/downstream/biomart.yaml'
    script:
        "../../scripts/figures/normTCreads_per-sample_plots.R"
