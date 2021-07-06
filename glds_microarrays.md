# glds_microarrays.R

## Description

This executable R script is the main function for processing GeneLab microarray datasets. It may be run with local data if runsheet and raw data files are supplied or called by glds_microarrays.sh which will pull data from the GeneLab repository API. It currently supports all Affymetrix platforms and Agilent one and two channel platforms. Support for Illumina Beadchip, Nimblegen one and two channel platforms as well as generic two channel platforms is in development.


## Required Packages

This software was written in R 4.0 and requires the following libraries to be installed:

Package | Source | Description
------- | ------ | -----------
optparse | Cran | Parses command line input options as an R object
rmarkdown | Cran | Convert R Markdown documents into a variety of formats
shiny | Cran | Builds interactive HTML documents
tidyverse | Cran | Collection of packages for tabular data manipulation and graphical output
BiocManager | Cran | Install Bioconductor sourced packages within an R session
oligo | Bioconductor | Pre-processing of Affymetrix and Nimblegen formatted raw array data
limma | Bioconductor | Pre-processing of Agilent, GenePix, and two-channel raw array data as well as linear model fitting and differential expression analysis
beadarray | Bioconductor | Pre-processing of Illumina Beadchip expression raw data
GEOquery | Bioconductor | Parsing Gene Expression Omnibus formatted metadata files

Cran packages can be installed from the R console like:

install.packages("BiocManager")

Bioconductor packages can be installed from the R console like:

BiocManager::install("oligo")

## Local Execution

From the command line run:
```bash
RSscript --vanilla <options>
```
Example:
```bash
Rscript --vanilla glds_microarrays.R --glds GLDS-205 --reports --runsheet ../data/GLDS-205/Metadata/GLDS-205_runsheet_based_on_lbl.csv ```
```

## API Execution

From the command line run:
```bash
sh glds_microarrays.sh <GLDS-XXX>
```

Example:
```bash
sh glds_microarrays.sh GLDS-205
```

### Input

Option | Type | Description
------ | ---- | -----------
--glds | character | Accession ID of the GLDS study to be processed
--reports | flag | If this flag is given, the script will generate raw and normalized quality assessment reports as HTML files
--runsheet | path | Path to GeneLab generated run_sheet.csv metadata file for local processing
--probe | path | Path to custom probe annotation file, primarily for two-channel datasets. If not supplied, the script will attempt to pull probe annotations from raw data files or from Bioconductor source databases
--out | path | Optional path to desired output directory. If not supplied, output will be copied to Processed_Data in the execution directory

Additional options are available for customized runs from local data utilizing ISAtab formatted metadata
Additional Options | Type | Description
------------------ | ---- | -----------
--files | path | Shell path to raw data file set. Accepts wildcards for multiple entries
--isa | path | Path to ISA.zip formatted metadata file
--species | character | Specify organism (e.g. "Arabidopsis thaliana")
--platform | character | Specify microarray platform (e.g. "Affymetrix")
--annotation | path | Path to organism gene level annotation file
--dir | path | Path to GeneLab curated directory structure top for a GLDS study

### Output

The script generates an output directory with the following structure:

00-RawData
01-NormalizedData
02-Limma_DGE
Metadata

00-RawData Files | Description
---------------- | -----------
raw files | Multiple raw data files with platform specific extensions (e.g. CEL, raw.txt, GPR, XYS, etc.)
raw_qa.html | Raw data quality assessment report generated with --reports flag
visualization_PCA_table.csv | PCA plotting table for data visualization of raw data
md5sum.txt | Checksum table for all raw data files

01-NormalizedData Files | Description
----------------------- | -----------
normalized.txt | Normalized expression values with probe identifiers
normalized-annotated.txt | Normalized expression values with gene level identifiers
normalized-annotated.rda | Normalized expression values with gene level identifiers in R object file format
normalized_qa.html | Normalized data quality assessment report generated with --reports flag
visualization_PCA_table.csv | PCA plotting table for data visualization of normalized data

02-Limma_DGE Files | Description
------------------ | -----------
contrasts.csv | Table of DGE group contrasts performed
differential_expression.csv | Gene level DGE results including annotations, expression values, group statistics, log fold change, p-value, and BH adjusted p-value for each group contrast
visualization_output_table.csv | Additional precalculated values to differential_expression.csv for rapid data visualization

Metadata Files | Description
-------------- | -----------
run_sheet.csv | GeneLab generated input file containing sample to file relasionships and other necessary metadata
probe_annotations.txt | Generalized probe annotation file constructed from source formatted annotation files
ISA.zip | ISAtab formatted metadata copied if using glds_microarrays.sh API script or if specified from a local source
  

## Platform Processing Functions

### affymetrix.R

This function takes the glds_microarray.R options object as input and generates output files customized to the Affymetrix platform. 
Preprocessing is performed with the Oligo package. RMA background correction and Quantiles normalization are performed to generate normalized expression values. For ST class arrays, core gene level estimation is performed. DGE analysis is performed with the Limma package. Control probes and unannotated probes are filtered. Multiple probes mapping to a gene identifier are aggregated by maximum interquartile range. Normalized expression values are fit to a linear model with empirical Bayes estimation of priors. Group contrast log fold change, t-test p-values, and Benjamini Hochberg adjusted p-values with FDR of 0.05 are calculated.

### agilent_one_channel.R

This function takes the glds_microarray.R options object as input and generates output files customized to the Affymetrix platform. 
Preprocessing is performed with the Oligo package. RMA background correction and Quantiles normalization are performed to generate normalized expression values. For ST class arrays, core gene level estimation is performed. DGE analysis is performed with the Limma package. Control probes and unannotated probes are filtered. Multiple probes mapping to a gene identifier are aggregated by maximum interquartile range. Normalized expression values are fit to a linear model with empirical Bayes estimation of priors. Group contrast log fold change, t-test p-values, and Benjamini Hochberg adjusted p-values with FDR of 0.05 are calculated.

### agilent_two_channel.R

This function takes the glds_microarray.R options object as input and generates output files customized to the Affymetrix platform. 
Preprocessing is performed with the Oligo package. RMA background correction and Quantiles normalization are performed to generate normalized expression values. For ST class arrays, core gene level estimation is performed. DGE analysis is performed with the Limma package. Control probes and unannotated probes are filtered. Multiple probes mapping to a gene identifier are aggregated by maximum interquartile range. Normalized expression values are fit to a linear model with empirical Bayes estimation of priors. Group contrast log fold change, t-test p-values, and Benjamini Hochberg adjusted p-values with FDR of 0.05 are calculated.

### illumina_expression.R

This function takes the glds_microarray.R options object as input and generates output files customized to the Affymetrix platform. 
Preprocessing is performed with the Oligo package. RMA background correction and Quantiles normalization are performed to generate normalized expression values. For ST class arrays, core gene level estimation is performed. DGE analysis is performed with the Limma package. Control probes and unannotated probes are filtered. Multiple probes mapping to a gene identifier are aggregated by maximum interquartile range. Normalized expression values are fit to a linear model with empirical Bayes estimation of priors. Group contrast log fold change, t-test p-values, and Benjamini Hochberg adjusted p-values with FDR of 0.05 are calculated.

## Quality Assessment Functions

### qa_summary_raw.rmd

This function takes platform specific raw data objects and generates an html report file containing chip level imageplots, intensity boxplots, density distribution plots, and PCA plots (for one channel datasets).

### qa_summary_normalized.rmd

This function takes platform specific normalized data objects and generates an html report file containing chip level intensity boxplots, density distribution plots, and PCA plots (for one channel datasets).
