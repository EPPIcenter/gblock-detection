// Process for getting spikein counts
process GET_SPIKEIN_COUNTS {
    tag "Getting Spikein Counts"
    publishDir(
        path: "${params.outDIR}",
        mode: 'copy'
    )
    label 'process_low'

    input:
    // The remaining reads to analyze; these fastqs should have
    // all panel reads removed at this point
    path (demultiplexed_spikeins_fastqs)

    output:
    path ('spikein_counts/final_spikein_counts.csv'), emit: spikein_counts

    script:
    """
    bash get-counts.sh \
        -c ${params.spikein_csv} \
        -d ${demultiplexed_spikeins_fastqs} \
        -o 'spikein_counts' \
        -t ${task.cpus}
    """
}