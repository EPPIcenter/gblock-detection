#!/usr/bin/env bash

cat ../results/stats/spikein_data.txt | tr "\t" "," > spikein_data.csv
cat ../results/amplicon_coverage.txt | tr "\t" "," > amplicon_coverage.csv
