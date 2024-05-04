#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Expand user directory if exists
outDIR = "${params.outDIR}".replaceFirst("^~", System.getProperty("user.home"))
readDIR = "${params.readDIR}".replaceFirst("^~", System.getProperty("user.home"))

// Boilerplate parameters
params.reads           = "${readDIR}/*_R{1,2}*.fastq.gz"
params.amplicon_info   = "$projectDir/resources/${params.target}/${params.target}_amplicon_info.tsv"
params.spikein_csv     = "$projectDir/resources/gblocks.csv"
params.primers_csv     = "$projectDir/resources/primers.tsv"

// Define parameters (can also be moved to nextflow.config)
include { DEMULTIPLEX_AMPLICONS } from './workflows/demultiplex_amplicons.nf'
include { SPIKEIN_ANALYSIS } from './workflows/spikein_analysis.nf'

workflow {
    // Create read pairs channel from fastq data
    read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )

    // Trim and demultiplex amplicons by amplicon
    DEMULTIPLEX_AMPLICONS(read_pairs)

    // Detect spikeins and create QC plots
    SPIKEIN_ANALYSIS(
        DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
        DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch,
        DEMULTIPLEX_AMPLICONS.out.demux_fastqs_ch,
        DEMULTIPLEX_AMPLICONS.out.unknown_fastqs_ch
    )
}
