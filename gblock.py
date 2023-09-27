import os
import argparse
import concurrent.futures
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm
import gzip
import re
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(message)s")

def generate_kmers(sequence, k):
    """Generates k-mers from a sequence."""
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]

def parse_gblock_sequences(filepath):
    gblocks = {}
    with open(filepath, 'r') as file:
        next(file)  # Skip header
        for line in file:
            line = line.strip().split(',')
            gblocks[line[0]] = line[1]
    return gblocks

def sanitize_filename(filename):
    """Extract the unique identifier from the filename."""
    # Use regex to extract the pattern
    match = re.search(r"Spike-In-(.+?)_R[12]_001.fastq.gz", filename)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Invalid filename pattern: {filename}")

def process_fastq(fastq_file, reads_directory, kmer_index, gblocks, kmer_size):
    filepath = os.path.join(reads_directory, fastq_file)
        
    counts = {name: 0 for name in gblocks.keys()}
    with gzip.open(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequence = str(record.seq)
            for kmer in generate_kmers(sequence, kmer_size):
                if kmer in kmer_index:
                    for name in kmer_index[kmer]:
                        counts[name] += 1
        
    sanitized_name = sanitize_filename(fastq_file)
    return sanitized_name, counts

def main(reads_directory, gblock_sequences, output_csv, kmer_size):
    logging.info(f"Using kmer-size: {kmer_size}")

    logging.info("Parsing GBlock sequences...")
    gblocks = parse_gblock_sequences(gblock_sequences)

    logging.info("Indexing GBlock k-mers...")
    kmer_index = {}
    for name, sequence in tqdm(gblocks.items(), desc="Indexing", ncols=100):
        for kmer in generate_kmers(sequence, kmer_size):
            if kmer not in kmer_index:
                kmer_index[kmer] = set()
            kmer_index[kmer].add(name)

    fastq_files = [f for f in os.listdir(reads_directory) if f.endswith('.fastq.gz')]
    logging.info(f"Processing {len(fastq_files)} fastq files...")

    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.cores) as executor:
        future_to_file = {executor.submit(process_fastq, fastq_file, reads_directory, kmer_index, gblocks, kmer_size): fastq_file for fastq_file in fastq_files}
        for future in tqdm(concurrent.futures.as_completed(future_to_file), total=len(fastq_files), desc="Processing fastq files", ncols=100):
            sanitized_name, counts = future.result()
            results[sanitized_name] = counts

    # Create DataFrame
    df = pd.DataFrame(results).fillna(0).astype(int).transpose()
    
    # Sort DataFrame alphanumerically by row index (sample names) and columns (gblock names)
    df = df.sort_index(axis=0, key=lambda x: x.str.lower()).sort_index(axis=1, key=lambda x: x.str.lower())
    
    # Write DataFrame to CSV
    print(df)
    df.to_csv(output_csv)

    # Display heatmap
    plt.figure(figsize=(25, len(fastq_files) * 0.6))
    sns.heatmap(np.log1p(df), cmap="coolwarm", cbar=True, xticklabels=True, yticklabels=True)
    plt.title(f"GBlock Sequence Hits in FASTQ Files (k-mer size: {kmer_size})")
    cbar_ax = plt.gcf().axes[-1]
    cbar_ax.set_title("log1p(counts)")
    plt.savefig("heatmap.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align fastq files to gblock sequences using k-mers and visualize hits.")
    parser.add_argument("--reads_directory", help="Directory containing fastq files.")
    parser.add_argument("--gblock_sequences", help="CSV file containing gblock sequences.")
    parser.add_argument("--output_csv", default="results.csv", help="Output CSV filename for the results.")
    parser.add_argument("--cores", type=int, default=4, help="Number of cores to use. Default is all available cores.")
    parser.add_argument("--kmer-size", type=int, default=125, help="Kmer size to use. This should be a value that allows you to identify GBlock sequences precisely with some tolerance for sequencing errors.")
    args = parser.parse_args()

    main(args.reads_directory, args.gblock_sequences, args.output_csv, args.kmer_size)
