/*
 * STEP - ALIGN-TO-REFERENCE
 * 
 */

process ALIGN_TO_REFERENCE {

  label 'process_high'

  input:
  path clusters
  path spikein_csv
  
  output:
  path("alignments.txt"), emit: alignments

  script:
  
  """
  Rscript ${projectDir}/bin/align_to_reference.R \
    --clusters ${clusters} \
    --refseq-csv ${spikein_csv} \
    --n-cores ${task.cpus}
  """
}