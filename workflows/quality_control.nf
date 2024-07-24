// workflows/quality_control.nf

include { PREPROCESS_COVERAGE } from '../modules/local/preprocess_coverage.nf'
include { QUALITY_REPORT } from '../modules/local/quality_report.nf'

workflow QUALITY_CONTROL {

    // Define inputs
    take:
    sample_coverage_files
    amplicon_coverage_files

    main:

    // Assuming 'sampleFiles' and 'ampliconFiles' are your channels for individual files
    sample_combined = sample_coverage_files.collect()
    amplicon_combined = amplicon_coverage_files.collect()

    // Initial Preprocessing
    PREPROCESS_COVERAGE(
        sample_combined,
        amplicon_combined
    )

    // If postprocessing coverage is provided, run the postprocessing workflow
    sample_coverage_ch = PREPROCESS_COVERAGE.out.sample_coverage
    amplicon_coverage_ch = PREPROCESS_COVERAGE.out.amplicon_coverage

    // Reporting
    QUALITY_REPORT(
        sample_coverage_ch,
        amplicon_coverage_ch,
        params.amplicon_info
    )
}
