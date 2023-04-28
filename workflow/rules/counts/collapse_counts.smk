def _input_bedFiles(wildcards):
    return [expand('resources/external/gencode_{release}/{genome}.transcripts.bed', release=GENCODE_RELEASE, genome=GENOME),
            expand('resources/external/gencode_{release}/{genome}.genes.bed', release=GENCODE_RELEASE, genome=GENOME)]


rule tx2genecounts:
    input: 
        unpack(_input_bedFiles),
        tcounTSV = rules.slam_count.output
    output:
        genecountsRData = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.RData',
        genecountsTSV = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.tsv',
        genecountsExtTSV = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.genecount.extended.tsv',
    resources:
        mem_mb = 4000
    log: 
        'logs/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_genecounts.log'
    conda:
        '../../envs/downstream/biomart.yaml'
    script:
        "../../scripts/counts/tx2genecounts.R"