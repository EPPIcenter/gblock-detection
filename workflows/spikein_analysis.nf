include { GET_SPIKEIN_COUNTS } from '../modules/local/get-spikein-counts.nf'
include { PLOT_SPIKEIN_COUNTS } from '../modules/local/plot-spikein-counts.nf'
include { SPIKEIN_TRIM } from '../modules/local/spikein-trim.nf'
include { CREATE_PRIMER_FILES } from '../modules/local/create_primer_files.nf'
include { ALIGN_TO_REFERENCE } from '../modules/local/align_to_reference.nf'
include { CREATE_SPIKEIN_TABLE } from '../modules/local/create_spikein_table.nf'

workflow SPIKEIN_ANALYSIS {

  take: 
  sample_summary_ch
  amplicon_summary_ch
  demux_fastqs_ch
  unknown_fastqs_ch

  main:

  CREATE_PRIMER_FILES(params.primers_csv)

  // demutltiplex spikeins
  SPIKEIN_TRIM(
    CREATE_PRIMER_FILES.out.fwd_primers,
    CREATE_PRIMER_FILES.out.rev_primers,
    unknown_fastqs_ch
  )

  // get the spikein counts from the trimmed spikin fastqs
  GET_SPIKEIN_COUNTS(
    SPIKEIN_TRIM.out.spikeins_demultiplexed
  )

  // plot the spikein counts in a heatmap
  PLOT_SPIKEIN_COUNTS(
    GET_SPIKEIN_COUNTS.out.spikein_counts.collect()
  )

  // create a spikein table
  CREATE_SPIKEIN_TABLE(
    SPIKEIN_TRIM.out.spikeins_demultiplexed.collect()
  )

  if (params.should_align_to_reference) {
    // align against our known spikins
    ALIGN_TO_REFERENCE(
      CREATE_SPIKEIN_TABLE.out.spikein_data,
      params.spikein_csv,
      
    )
  }
}