# Use an official R base image from the Rocker project
FROM rocker/r-ubuntu:22.04

# Install necessary packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-pip \
    libssl-dev \
    libbz2-dev \
    libsdl1.2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    bwa \
    bcftools \
    samtools \
    pandoc \
    jq \
    flash \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install cutadapt==4.4 pandas==2.0.3 bio==1.5.9

# Set R configurations
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/ && \
    echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | \
    tee /usr/local/lib/R/etc/Rprofile.site | \
    tee /usr/lib/R/etc/Rprofile.site

# Install R packages using Rscript
RUN Rscript -e 'install.packages("remotes", version = "2.5.0")' && \
    Rscript -e 'remotes::install_cran(c("tidyverse", "reshape2", "argparse", "ggbeeswarm", "gridExtra", "rmarkdown", "foreach", "doMC"), upgrade="never", version=c("2.0.0", "1.4.4", "2.1.6", "0.6.1", "2.3", "2.17", "1.5.2", "1.3.8"))'

# Install Bioconductor Packages
RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager"); }; BiocManager::install(version = "3.18");'

RUN Rscript -e 'BiocManager::install("dada2", version = "3.18", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("muscle", version = "3.18", ask = FALSE)'
RUN Rscript -e 'BiocManager::install("BSgenome", version = "3.18", ask = FALSE)'

# Define default command to start bash
CMD ["/bin/bash"]
