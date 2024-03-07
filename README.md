# NRCD_analysis

Code and data form the "Microbiome gut community structure and functionality is associated with symptoms severity in non-responsive celiac disease patients undergoing gluten-free diet" paper.

Raw data can be found in ENA under accession number PRJEB65879. [Here](https://www.biorxiv.org/content/10.1101/2023.10.14.562350v2) is a preprint of the study published in bioRxiv.

## Table of contents
* [MetOrigins](#data)
* [Scripts](#scripts)

## MetOrigins
This directory contains the results after performing [MetOrigin](https://metorigin.met-bioinformatics.cn/home/) analysis.

## R Scripts
This directory contains the scripts used for upstream analysis in R.
- `phyloseq_analysis.R`: script used to create `phyloseq` object from metaphlan output.
- `MFA.R`: script used to perform MFA analysis and select informative variables. 
- `Networks.R`: script used to construct co-occurrence networks. 
- `complex_heatmap.R`: script used to perform Pearson correlation analysis and visualization.
