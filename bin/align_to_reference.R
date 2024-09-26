library(argparse)

parser <- ArgumentParser(description='Aligns sequences (DADA2 clusters) against a reference of known sequences')
parser$add_argument('--clusters', type="character", required=TRUE, help="RDS Clusters from DADA2. This is the main output from the DADA module.")
parser$add_argument('--refseq-fasta', type="character")
parser$add_argument('--alignment-threshold', type="integer", default = 60)
parser$add_argument('--n-cores', type = 'integer', default = 1, help = "Number of cores to use. Ignored if running parallel flag is unset.")
parser$add_argument('--amplicon-table', type="character", required=TRUE, help = "Amplicon table with primer pools. This is used to organize the sequence table by amplicon.")

args <- parser$parse_args()
print(args)

library(stringr)
library(dplyr)
library(dada2)
library(foreach)
library(parallel)
library(BSgenome)
library(tidyr)
library(doMC)
library(tibble)
library(ggplot2)
library(Biostrings)
library(magrittr)

## FOR DEUBGING
# setwd("/home/eppicenter/Desktop")
# args=list()
# args$clusters="spikeins_data.txt"
# args$refseq_csv="spikeins.csv"
# args$parallel=TRUE
# args$n_cores=2

clusters=read.table(args$clusters, header=T)

## DELETE BELOW
spikin_csv <- read.csv(args$refseq_csv)
spikins <- spikin_csv$GBlockSequence
names(spikins) <- spikin_csv$SpikeinID
ref_sequences <- DNAStringSet(spikins)

# register number of cores to use
registerDoMC(args$n_cores)

# DADA2 should have corrected any substitution errors during sequencing. Therefore, any 
# mismatches are expected to be real. A small penalty is added when there is a mismatching base
# to catch large differences that are likely not real. Gap penalties are also set below to 
# help filter sequences that are truly different from the reference. 
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)

# A06_0SI_1:10PR_042024_SI_IX
# B06_0SI_1:10PR_042024_SI_IX
# C06_0SI_1:10PR_042024_SI_IX
# D06_0SI_1:10PR_042024_SI_IX
# E06_0SI_1:10PR_042024_SI_IX
# F06_0SI_1:10PR_042024_SI_IX

#########  IMPORTANT: NON PF SEQUENCES NEED TO BE INCLUDED IN REFERENCE!
refseq.seq1 <- spikin_csv %>%
  filter(SpikeinID == )

# This object contains the aligned ASV sequences
df_aln <- foreach(seq1 = 1:nrow(clusters), .combine = "bind_rows") %dopar% {
  # the alignment is performed only with the reference sequence for the corresponding locus as we have this information from the demultiplexing step
  aln <- pairwiseAlignment(refseq.seq1, str_remove_all(clusters$Spikein[seq1],"N"), substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
  patt <- c(alignedPattern(aln), alignedSubject(aln))
  ind <- sum(str_count(as.character(patt),"-"))
  data.frame(
    SampleID = clusters$SampleID[seq1],
    asv = clusters$asv[seq1],
    hapseq = as.character(patt)[2],
    refseq = as.character(patt)[1],
    refid = clusters$locus[seq1],
    score = score(aln),
    indels = ind
  )
}


write.table(df_aln,file="alignments.txt",quote=F,sep="\t",col.names=T,row.names=F)
