# Bioinformatics pipeline for Agilent two-channel microarray platform

> **This document holds an overview and some example commands of how GeneLab processes Agilent microarray datasets. Exact processing commands for specific datasets that have been released is available in this repository [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory and is also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** July 7, 2021  
**Revision:** -  
**Document Number:**   

**Submitted by:**  
Homer Fogle (GeneLab Analysis Team)  

**Approved by:**  


---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**Pre-processing**](#pre-processing)
    - [1. Raw Data QA](#1-raw-data-qa)
    - [2. Background correction and normalization](#2-background-correction)
    - [3. Normalized Data QA](#3-normalized-data-qa)
    - [4. Probe filtering](#4-probe-filtering)
    - [5. Probe annotation](#5-probe-annotation)
  - [**Differential Gene Expression**](#differential-gene-expression)
    - [6. Construct Model Matrix](#6-model-matrix)
    - [7. Fit Linear Model](#7-linear-model)
    - [8. Construct Contrast Matrix](#8-contrast-matrix)
    - [9. Fit Contrast Model](#9-contrast-model)
    - [10. Perform Empirical Bayes Moderation](#10-empirical-bayes)
    - [11. Generate DGE Statistics](#11-dge-stats)


---

# Software used

|Program|Version*|Relevant Links|
|:------|:-----:|------:|
|R|4.1|[https://cran.r-project.org/](https://cran.r-project.org/)|
|optparse|1.6.6|[https://cran.r-project.org/web/packages/optparse/index.html](https://cran.r-project.org/web/packages/optparse/index.html)|
|BiocManager|1.30|[https://cran.r-project.org/web/packages/BiocManager/index.html](https://cran.r-project.org/web/packages/BiocManager/index.html)|
|tidyverse|1.3.1|[https://cran.r-project.org/web/packages/tidyverse/index.html](https://cran.r-project.org/web/packages/tidyverse/index.html)|
|rmarkdown|2.9|[https://cran.r-project.org/web/packages/rmarkdown/index.html](https://cran.r-project.org/web/packages/rmarkdown/index.html)|
|shiny|1.6|[https://cran.r-project.org/web/packages/shiny/index.html](https://cran.r-project.org/web/packages/shiny/index.html)|
|oligo|3.13|[https://bioconductor.org/packages/release/bioc/html/oligo.html](https://bioconductor.org/packages/release/bioc/html/oligo.html)|
|limma|3.48|[https://bioconductor.org/packages/release/bioc/html/limma.html](https://bioconductor.org/packages/release/bioc/html/limma.html)|
|GEOquerry|2.60|[https://bioconductor.org/packages/release/bioc/html/GEOquery.html](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)|
|beadarray|2.42|[https://bioconductor.org/packages/release/bioc/html/beadarray.html](https://bioconductor.org/packages/release/bioc/html/beadarray.html)|
|qa_summary_raw.rmd|1.0|[https://bioconductor.org/packages/release/bioc/html/GEOquery.html](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)|
|qa_summary_normalized.rmd|1.0|[https://bioconductor.org/packages/release/bioc/html/beadarray.html](https://bioconductor.org/packages/release/bioc/html/beadarray.html)|
|genefilter|1.74|[https://bioconductor.org/packages/release/bioc/html/genefilter.html](https://bioconductor.org/packages/release/bioc/html/genefilter.html)|
|annotationDbi|1.54|[https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)|
>**\*** Exact versions utilized for a given dataset are available along with the processing commands for each specific dataset (this is due to how the system may need to be updated regularly).

---

# General processing overview with example commands

> Exact processing commands for specific datasets is available in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory of this repository, as well as being provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

## Pre-processing

### 1. Raw Data QA

```
rmarkdown::render("qa_summary_raw.Rmd","html_document", output_file="raw_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"00-RawData"))

```

**Parameter Definitions:**

* `input` – input RMD file to render
* `output_format` – specify rendering output file format (e.g. HTML, PDF)
* `output_file` - specify filename of output file
* `output_dir` - specify directory of output file

**Input data:**

* opt$glds (R character object containing GLDS study identifier)
* raw (R oligoclass object containing raw data)

**Output data:**

* *qa_raw.html (QA output html summary of raw data)



#### 2. Background Correction and Normalization

```
data <- oligo::rma(raw, target = "core", background=TRUE, normalize=TRUE)
```

**Parameter Definitions:**

*	`object` – R oligoclass object containing raw data
*	`target` – Level of summarization (only for Exon/Gene arrays)
*	`background` – Logical - perform RMA background correction?
*	`normalize` – Logical - perform quantile normalization?

**Input data:**

* raw (R oligoclass object containing raw data)

**Output data:**

* data (R oligoclass object of normalized data)


<br>  

---

### 3. Normalized Data QA

```      
rmarkdown::render("qa_summary_normalized.Rmd","html_document", output_file="normalized_qa",output_dir=file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))

```

**Parameter Definitions:**

* `input` – input RMD file to render
* `output_format` – specify rendering output file format (e.g. HTML, PDF)
* `output_file` - specify filename of output file
* `output_dir` - specify directory of output file

**Input data:**

* opt$glds (R character object containing GLDS study identifier)
* data (R oligoclass object containing normalized data)

**Output data:**

* *qa_raw.html (QA output html summary of raw data)

<br>

---

### 4.0 Filter Probes
```
data.filt <- genefilter::nsFilter(data, require.entrez=TRUE,
#                                        remove.dupEntrez=TRUE, var.func=IQR,
#                                        var.cutoff=0.5, var.filter=TRUE,
#                                        filterByQuantile=TRUE, feature.exclude="^AFFX")
```

**Parameter Definitions:**

*	`eset` – an ExpressionSet object 
*	`require.entrez` – 	If TRUE, filter out features without an Entrez Gene ID annotation.
*	`remove.dupEntrez` - If TRUE and there are features mapping to the same Entrez Gene ID (or equivalent), then the feature with the largest value of var.func will be retained and the other(s) removed.
*	`var.func` - The function used as the per-feature filtering statistic.
*	`var.cutoff` - A numeric value. If var.filter is TRUE, features whose value of var.func is less than either: the var.cutoff-quantile of all var.func values (if filterByQuantile is TRUE), or var.cutoff (if filterByQuantile is FALSE) will be removed.
*	`var.filter` - A logical indicating whether to perform filtering based on var.func.
*	`filterByQuantile` - 	A logical indicating whether var.cutoff is to be interprested as a quantile of all var.func values (the default), or as an absolute value.
* `feature.exclude` - A character vector of regular expressions. Feature identifiers (i.e. value of featureNames(eset)) that match one of the specified patterns will be filtered out. The default value is intended to filter out Affymetrix quality control probe sets.


**Input data:**

* data (R oligoclass ExpressionSet object)

**Output data:**

* data.filt (R oligoclass ExpressionSet object)



#### 5. Probe annotation
```
annotation<-mapIds(database,keys = keys,keytype = keytype, column = column,multiVals = "first")
```

**Parameter Definitions:**

*`database` – the AnnotationDb object
*`keys` – the keys to select records for from the database. 
*`keytype` – the keytype that matches the keys used. For the select methods, this is used to indicate the kind of ID being used with the keys argument
*`column` – the column to search on (for mapIds). 
*`multiVals` - This value means that when there are multiple matches only the 1st thing that comes back will be returned. 

**Input data:**

* database (R AnnotationDbi object database specific to the Affynetrix platforn)
* keys (R dataframe of probe IDs pulled from raw data object
* keytype (character string specifying database keytype
* column (character string specifying database identifier to pull

**Output data:**

* annotation (R dataframe of selected database identifiers)


<br>

---

## Differential Gene Expression
### 6. Construct Model Matrix
```
design <- model.matrix(formula)
```

**Parameter Definitions:**  

*`formula` – R formula object

**Input data:**

* ~0 + groups (R formula with groups dataframe ordered by array)

**Output data:**

* design (design matrix from the description given in terms(formula object))


<br>

---

### 7. Fit Linear Model

```
fit <- limma::lmFit(object, design)
```

**Parameter Definitions:**  

*	`object` – A matrix-like data object containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples
*	`design` – the design matrix of the microarray experiment


**Input data:**

* data.filt (filtered normalized expression data)
* design (model matrix of experiment

**Output files:**

* fit (an MArrayLM object)


### 8. Construct contrast matrix

```
cont.matrix <- makeContrasts(contrasts = contrasts,levels=design)
```

**Parameter Definitions:**  

*`contrasts` – 	character vector specifying contrasts

*`levels`– design matrix


**Input data:**

* contrasts (group contrast character vector in R formula format)
* design (design matrix object

**Output files:**

* contr.matrix (matrix of contrast conditions)

<br>

---

### 9 Fit Contrast model
```
contrast.fit <- contrasts.fit(fit, cont.matrix)
```
**Parameter Definitions:**

*	`fit` – an MArrayLM object or a list object produced by the function lm.series or equivalent

*	`contrasts` – numeric matrix with rows corresponding to coefficients in fit and columns containing contrasts.


**Input data:**

* fit (linear model fit object)
* cont.matrix (contrasts matrix)

**Output data:**

* contrast.fit (contrast fit model object)


<br>

---

### 10. Empirical Bayes Moderation

```
contrast.fit <- eBayes(contrast.fit)
```

**Parameter Definitions:**
*	`fit` – an MArrayLM fitted model object produced by lmFit or contrasts.fit.


**Input data:**

* contrast.fit (contrast model fit object)


**Output data:**

* contrast.fit (contrast model fit object)


#### 11. Generate DGE Statistics
```
top <- topTable(contrast.fit, coef = i, number = Inf, genelist = contrast.fit$genes$ID, adjust.method = "BH", sort.by = "none")
```

**Parameter Definitions:**  

* `fit` – an object of class MArrayLM as produced by lmFit and eBayes.
* `coef` – column number or column name specifying which coefficient or contrast of the linear model is of interest.
* `number` - maximum number of genes to list
* `genelist`- data frame or character vector containing gene information.
* `adjust.method` - method used to adjust the p-values for multiple testing.
* `sort.by` - character string specifying which statistic to rank the genes by.


**Input data:**

* contrast.fit (linear model object of group contrasts)

**Output data:**

* top (R dataframe of gene level DGE statistics)

<br>

---


---
