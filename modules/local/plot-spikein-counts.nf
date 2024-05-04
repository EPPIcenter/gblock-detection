// Process for generating plots
process PLOT_SPIKEIN_COUNTS {
    tag "Generating Plots"
    publishDir(
        path: "${params.outDIR}/plots",
        mode: 'copy'
    )
    label 'process_single'

    input:
    path (counts_files, stageAs: "counts_files?")

    output:
    path 'spikein_counts_heatmap.png', emit: spikein_counts_heatmap

    script:
    """
    Rscript ${projectDir}/bin/spikein_detection_heatmap.R \
        --input ${counts_files} \
        --output spikein_counts_heatmap.png
    """
}