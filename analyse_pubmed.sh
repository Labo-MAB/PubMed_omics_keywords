#!/bin/bash

# This script is only a wrapper around 2 adapted tools to analyse the pubmed
# corpus data. The result file generated by pubmed_trend_analysis can then be
# analysed with create_heatmaps.r and create_wordclouds.py to generate the plots

# Orginal version of the tools can be found at:
# https://github.com/lab42open-team/pubmed_trend_analysis
# https://github.com/sahansera/medline-pubmed-extractor

# Download and unzip the pubmed data
mkdir -p pubmed/xml/
for i in $(seq -w 0001 1219)
do
    wget https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed24n${i}.xml.gz -O pubmed/xml/pubmed24n${i}.xml.gz
    gzip -d pubmed/xml/pubmed24n${i}.xml.gz
done

# Convert xml files to tsv files using medline-pubmed-extractor
mkdir -p pubmed/tsv/
tools/medline-pubmed-extractor/MedlineExtractor/bin/Release/MedlineExtractor pubmed/xml/ pubmed/tsv/
gzip pubmed/tsv/*

# Extract the articles containing the keywords using pubmed_trend_analysis
tools/pubmed_trend_analysis/dig_analysis.sh -k omics_keywords.txt -d pubmed/tsv/  -p omics
