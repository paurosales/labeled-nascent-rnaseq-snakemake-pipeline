def _input_filter(wildcards):
    return TARGETS['slam_filter']

def _input_tcounts(wildcards):
    return TARGETS['slam_count']


rule alley_summary:
    input: 
        filteredBAM = unpack(_input_filter),
        counts = unpack(_input_summary)
    output:
        'results/quality_control/summary/slamseq_summary.txt'
    params:
        counts_dir = 'results/counts/'
    resources:
        mem_mb = 12000
    threads: 16
    conda:
        '../../envs/raw_processing/slamdunk.yaml'
    shell:
        """
            alleyoop summary -o {output} -t {params.counts_dir} \
            {input.filteredBAM}
        """