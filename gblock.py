import os
import argparse
import concurrent.futures
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm
import gzip
import re
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(message)s")

def align(query, gblocks):
    best_score = 0
    best_id = None
    best_errors = 0
    
    # Set up the alignment parameters
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    
    for ref_id, ref_seq in gblocks.items():
        alignments = aligner.align(query, ref_seq)
        # Check if there are any alignments
        if alignments:
            top_alignment = alignments[0]
            score = top_alignment.score
            aligned_query = top_alignment.aligned[0][0]
            aligned_ref = top_alignment.aligned[1][0]
            errors = sum(1 for a, b in zip(aligned_query, aligned_ref) if a != b)
            if score > best_score:
                best_score = score
                best_id = ref_id
                best_errors = errors

    return best_id, best_score, best_errors

def generate_kmers(sequence, kmer_size):
    """Generate k-mers from the sequence."""
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def create_kmer_index(gblocks, kmer_size):
    kmer_index = {}
    for gblock_name, gblock_sequence in gblocks.items():
        for kmer in generate_kmers(gblock_sequence, kmer_size):
            if kmer not in kmer_index:
                kmer_index[kmer] = []
            kmer_index[kmer].append(gblock_name)
    return kmer_index

def parse_gblock_sequences(filepath):
    gblocks = {}
    with open(filepath, 'r') as file:
        next(file)  # Skip header
        for line in file:
            line = line.strip().split(',')
            gblocks[line[0]] = line[1]
    return gblocks

def sanitize_filename(filename):
    match = re.search(r"(?:Spike-In|Control\d+)-(.+?)_R[12]_.*$", filename)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Invalid filename pattern: {filename}")

def process_fastq(fastq_file, reads_directory, kmer_index, gblocks, kmer_size, method="alignment"):
    filepath = os.path.join(reads_directory, fastq_file)
        
    counts = {name: 0 for name in gblocks.keys()}
    with gzip.open(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequence = str(record.seq)
            
            if method == "kmer":
                for kmer in generate_kmers(sequence, kmer_size):
                    if kmer in kmer_index:
                        for gblock_name in kmer_index[kmer]:
                            counts[gblock_name] += 1
            else:  # alignment
                best_id, _, _ = align(sequence, gblocks)
                if best_id:
                    counts[best_id] += 1
        
    sanitized_name = sanitize_filename(fastq_file)
    return sanitized_name, counts


def main(reads_directory, gblock_sequences, output_csv, kmer_size, cores, method):
    logging.info("Parsing GBlock sequences...")
    gblocks = parse_gblock_sequences(gblock_sequences)

    # Create the kmer_index if the method is set to kmer
    kmer_index = {}
    if method == "kmer":
        kmer_index = create_kmer_index(gblocks, kmer_size)

    fastq_files = [f for f in os.listdir(reads_directory) if f.endswith('.fastq.gz')]
    logging.info(f"Processing {len(fastq_files)} fastq files...")

    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        future_to_file = {executor.submit(process_fastq, fastq_file, reads_directory, kmer_index, gblocks, kmer_size, method): fastq_file for fastq_file in fastq_files}
        
        progress_bar = tqdm(total=len(fastq_files), desc="Processing fastq files", ncols=100)
        for future in concurrent.futures.as_completed(future_to_file):
            sanitized_name, counts = future.result()
            results[sanitized_name] = counts
            progress_bar.set_description(f"Retrieved results for: {sanitized_name}")
            progress_bar.update(1)
        progress_bar.close()

    df = pd.DataFrame(results).fillna(0).astype(int).transpose()
    df = df.sort_index(axis=0, key=lambda x: x.str.lower()).sort_index(axis=1, key=lambda x: x.str.lower())
    
    print(df)
    df.to_csv(output_csv)

    plt.rcParams['axes.labelsize'] = 26  
    plt.rcParams['xtick.labelsize'] = 24  
    plt.rcParams['ytick.labelsize'] = 24  
    plt.rcParams['axes.titlesize'] = 46   

    plt.figure(figsize=(25, len(fastq_files) * 0.6))
    sns.heatmap(np.log1p(df), cmap="coolwarm", cbar=True, xticklabels=True, yticklabels=True)
    
    # Customizing title based on method and kmer size
    if method == "kmer":
        title = f"GBlock Sequence Hits in FASTQ Files using K-mer method (k={kmer_size})"
    else:
        title = f"GBlock Sequence Hits in FASTQ Files using Alignment method"

    plt.title(title, fontsize=20)
    
    cbar_ax = plt.gcf().axes[-1]
    cbar_ax.set_title("log1p(counts)", fontsize=16)
    plt.tight_layout()  
    plt.savefig("heatmap.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align fastq files to gblock sequences and visualize hits.")
    parser.add_argument("--reads-directory", help="Directory containing fastq files.")
    parser.add_argument("--gblock-sequences", help="CSV file containing gblock sequences.")
    parser.add_argument("--output-csv", default="results.csv", help="Output CSV filename for the results.")
    parser.add_argument("--cores", type=int, default=4, help="Number of cores to use.")
    parser.add_argument("--kmer-size", type=int, default=100)
    parser.add_argument("--method", choices=["kmer", "alignment"], default="alignment", help="Method to use for matching sequences. Choices are: kmer, alignment. Default is alignment.")
    args = parser.parse_args()

    main(args.reads_directory, args.gblock_sequences, args.output_csv, args.kmer_size, args.cores, args.method)
