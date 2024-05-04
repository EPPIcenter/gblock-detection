# README

## To run

Remove adapters and demultiplex fastq files.

```bash
bash trim.sh -f /path/to/fastqs.gz -t threads -r results
```

Get spikein counts for each sample.

```bash
bash get-counts.sh
```

Get plots.

```bash
Rscript plot.R
```


## To use Snakemake

Install using:

```bash
pip install snakemake
```

This will run the `Snakefile` in the current directory. To run with additional cores, use:

```bash
snakemake --cores [number_of_cores]
```

### Use Snakemake with Apptainer

I have an apptainer image put together, which may be more streamlined for you rather than installing the dependencies yourself.

First make sure you have [apptainer]. Then run:

```bash
apptainer build spikein-analysis.sif Apptainer
```

