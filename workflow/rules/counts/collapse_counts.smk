def _input_ensembl_tabs(wildcards):
    return [expand('resources/external/ensembl/{genome}.geneset.tsv', genome=GENOME),
            expand('resources/external/ensembl/{genome}.transcriptset.tsv', genome=GENOME)]


rule gene2transcript_counts:
    input: 
        unpack(_input_ensembl_tabs),
        tcounTSV = rules.slam_count.output
    output:
        geneCountsRData = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.RData',
        geneCountsTSV = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.tsv',
        geneCountsExtraTSV = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.extended.tsv'
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
#         genesetTSV = rules.get_ensembl_geneset.output.genesetTSV,
#         transcriptsetTSV = rules.get_ensembl_geneset.output.transcriptsetTSV
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