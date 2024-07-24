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

include { QUALITY_CONTROL } from './workflows/quality_control.nf'

workflow {
    // Create read pairs channel from fastq data
    read_pairs = channel.fromFilePairs( params.reads, checkIfExists: true )

    // Trim and demultiplex amplicons by amplicon
    DEMULTIPLEX_AMPLICONS(read_pairs)

    // Create the quality report now
    QUALITY_CONTROL(
      DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
      DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch
    )

    // Detect spikeins and create QC plots
    SPIKEIN_ANALYSIS(
        DEMULTIPLEX_AMPLICONS.out.sample_summary_ch,
        DEMULTIPLEX_AMPLICONS.out.amplicon_summary_ch,
        DEMULTIPLEX_AMPLICONS.out.demux_fastqs_ch,
        DEMULTIPLEX_AMPLICONS.out.unknown_fastqs_ch
    )
}


workflow.onComplete {
  def outputDir = new File("${params.outDIR}/run")
  if (!outputDir.exists()) {
      outputDir.mkdirs()
  }

  record_params()
  record_runtime()
}


/* Record parmameters
 *
 * Records set parameters for the run
 *
 */
def record_params() {
    
    def output = new File("${params.outDIR}/run/parameters.tsv")

    params.each{ k, v -> 
        output.append("${k}\t${v}\n")
    }
}


/* Record runtime information
 *
 * Records runtime and environment information and writes summary to a tabulated (tsv) file
 *
 */
def record_runtime() {

    def output = new File("${params.outDIR}/run/runtime.tsv")

    // Append ContainerEngine first as shown in your example
    output.append("PipelineVersion\t${workflow.manifest.version}\n")
    output.append("ContainerEngine\t${workflow.containerEngine}\n")
    output.append("Duration\t${workflow.duration}\n")
    output.append("CommandLine\t${workflow.commandLine}\n")
    output.append("CommitId\t${workflow.commitId}\n")
    output.append("Complete\t${workflow.complete}\n")
    output.append("ConfigFiles\t${workflow.configFiles.join(', ')}\n")
    output.append("Container\t${workflow.container}\n")
    output.append("ErrorMessage\t${workflow.errorMessage}\n")
    output.append("ErrorReport\t${workflow.errorReport}\n")
    output.append("ExitStatus\t${workflow.exitStatus}\n")
    output.append("HomeDir\t${workflow.homeDir}\n")
    output.append("LaunchDir\t${workflow.launchDir}\n")
    output.append("Manifest\t${workflow.manifest}\n")
    output.append("Profile\t${workflow.profile}\n")
    output.append("ProjectDir\t${workflow.projectDir}\n")
    output.append("Repository\t${workflow.repository}\n")
    output.append("Resume\t${workflow.resume}\n")
    output.append("Revision\t${workflow.revision}\n")
    output.append("RunName\t${workflow.runName}\n")
    output.append("ScriptFile\t${workflow.scriptFile}\n")
    output.append("ScriptId\t${workflow.scriptId}\n")
    output.append("ScriptName\t${workflow.scriptName}\n")
    output.append("SessionId\t${workflow.sessionId}\n")
    output.append("Start\t${workflow.start}\n")
    output.append("StubRun\t${workflow.stubRun}\n")
    output.append("Success\t${workflow.success}\n")
    output.append("UserName\t${workflow.userName}\n")
    output.append("WorkDir\t${workflow.workDir}\n")
    output.append("NextflowBuild\t${nextflow.build}\n")
    output.append("NextflowTimestamp\t${nextflow.timestamp}\n")
    output.append("NextflowVersion\t${nextflow.version}\n")
}