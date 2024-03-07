# NRCD_analysis
Code and data form the "Microbiome gut community structure and functionality is associated with symptoms severity in non-responsive celiac disease patients undergoing gluten-free diet" paper.

Raw data can be found in ENA under accession number PRJEB65879. [Here](https://www.biorxiv.org/content/10.1101/2023.10.14.562350v2) is a preprint of the study published in bioRxiv.

## Table of contents
* [MetOrigins](#MetOrigins)
* [Metadata](#Metadata)
* [Scripts](#RScripts)

## MetOrigins
This directory contains the results after performing [MetOrigin](https://metorigin.met-bioinformatics.cn/home/) analysis.

## Metadata
This directory contains metadata files used for analysis.
- `Antrophometrics`: Include anthropometric data of the patients.
- `Biochemistry`: Includes biochemistry, haematology, inflammatory and intestinal mucosal damage markers.
- `CDAT`: Score of Celiac Dietary Adherence Test.
- `CD-QOL`: Celiac Disease Quality of Life (CD-QOL) survey results.
- `CeD-Pro`: Results of the Celiac Disease Patient Reported Outcome.
- `GSRT-CeD`: Results of the Gastrointestinal Symptoms Rating Scale. 

## R Scripts
This directory contains the scripts used for upstream analysis in R.
- `phyloseq_analysis.R`: script used to create `phyloseq` object from metaphlan output.
- `MFA.R`: script used to perform MFA analysis and select informative variables. 
- `Networks.R`: script used to construct co-occurrence networks. 
- `complex_heatmap.R`: script used to perform Pearson correlation analysis and visualization.
