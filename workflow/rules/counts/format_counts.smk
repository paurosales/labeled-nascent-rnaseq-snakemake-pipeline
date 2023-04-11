rule gene2transcript_counts:
    input: 
        tcounTSV = rules.slam_count.output,
        genesetTSV = rules.get_ensembl.output.genesetTSV,
        transcriptsetTSV = rules.get_ensembl.output.transcriptsetTSV
    output:
        geneCountsRData = 'results/{proj}/slamdunk/count/{label}_geneCounts.RData',
        geneCountsTSV = 'results/{proj}/slamdunk/count/{label}_tcount_by_gene.tsv',
        geneCountsExtraTSV = 'results/{proj}/slamdunk/count/{label}_tcount_by_gene_extended.tsv'
    params:
        ensembl_version = config['ensembl_version']
    resources:
        mem_mb = 4000
    log: 
        'logs/{proj}/slamdunk/{label}_collapse_counts.log'
    conda:
        '../../envs/biomart.yaml'
    script:
        "../../scripts/counts/collapse_transcripts.R"


# rule rates_collapse_transcripts:
#     input: 
#         convRatesTSV = rules.alley_utrrates.output,
#         genesetTSV = rules.get_ensembl.output.genesetTSV,
#         transcriptsetTSV = rules.get_ensembl.output.transcriptsetTSV
#     output:
#         convRatesGeneTSV = 'results/{proj}/alleyoop/utrrates/{label}_convrates_by_gene.tsv'
#     params:
#         genome_build = config['genome_build']
#     resources:
#         mem_mb = 4000
#     log: 
#         'results/{proj}/logs/slamdunk/{label}_collapse_rates.log'
#     conda:
#         '../envs/deseq2.yaml'
#     script:
#         "../scripts/rates_collapse_transcripts.R"