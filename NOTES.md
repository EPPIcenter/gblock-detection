# NOTES

## Feb 6 2024

Use the following command if cutadapt throws a "too many open files" error.

```bash
ulimit -S -n 4096
```

## May 3 2024

Use the following command to make the gblock fasta for a reference (and save it in resources)

```bash
awk -F',' 'NR>1 {print ">"$1"\n"$2}' spikein_data.csv > spikein_sequences.fasta
```