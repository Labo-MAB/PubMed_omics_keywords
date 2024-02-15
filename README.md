# Pubmed omic keywords analysis

 __Authors__: Francis Bourassa ([@francis-B](https://github.com/Francis-B))

 __Email__: <Francis.Bourassa@USherbrooke.ca>

## Description

This repository contains the scripts used to produce heatmaps and wordclouds of <u>ARTICLENAME</u>.

The PubMed Corpus files were converted to tsv files using a modified version of the [medline-pubmed-extractor](https://github.com/sahansera/medline-pubmed-extractor) tool developped by Sahan Serasinghe ([@sahansera](https://github.com/sahansera)) and the keywords were searched in these tsv with a modified version of the [pubmed_trend_analysis](https://github.com/lab42open-team/pubmed_trend_analysis) tool developped by [lab42OPEN](https://github.com/lab42open-team). Both modified versions can be found in *tool/*.

## How to run

To reproduce the analysis, one can use the *analyse_pubmed.sh* which download the Pubmed Corpus files and call the tools mentionned above. Finally, the output of pubmed_trend_analysis, which can be found in *results/*, can be analysed with the *create_heatmaps.R* and *create_wordclouds.py* scripts to reproduce the figures.

To ease the last step, the conda environnements used to run the R and python scripts can be found in *conda/*.
