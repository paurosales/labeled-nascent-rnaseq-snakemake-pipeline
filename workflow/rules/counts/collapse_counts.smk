rule gene2transcript_counts:
    input: 
        tcounTSV = rules.slam_count.output,
        genesetTSV = rules.get_ensembl.output.genesetTSV,
        transcriptsetTSV = rules.get_ensembl.output.transcriptsetTSV
    output:
        geneCountsRData = 'results/slamdunk/count/{sample_type}_{treatment}_Bio-rep_{bio_rep}_geneCounts.RData',
        geneCountsTSV = 'results/slamdunk/count/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.tsv',
        geneCountsExtraTSV = 'results/slamdunk/count/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.extended.tsv'
    params:
        ensembl_version = config['ENSEMBL']['VERSION']
    resources:
        mem_mb = 4000
    log: 
        'logs/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}_collapse_counts.log'
    conda:
        '../../envs/downstream/biomart.yaml'
    script:
        "../../scripts/counts/collapse_transcripts.R"


# rule rates_collapse_transcripts:
#     input: 
#         convRatesTSV = rules.alley_utrrates.output,
#         genesetTSV = rules.get_ensembl.output.genesetTSV,
#         transcriptsetTSV = rules.get_ensembl.output.transcriptsetTSV
#     output:
#         convRatesGeneTSV = 'results/alleyoop/utrrates/{sample_type}_{treatment}_Bio-rep_{bio_rep}_convrates_by_gene.tsv'
#     params:
#         genome_build = config['genome_build']
#     resources:
#         mem_mb = 4000
#     log: 
#         'results/logs/slamdunk/{sample_type}_{treatment}_Bio-rep_{bio_rep}_collapse_rates.log'
#     conda:
#         '../envs/deseq2.yaml'
#     script:
#         "../scripts/rates_collapse_transcripts.R"