# Bioinformatics pipeline for Affymetrix microarray platform

> **This document holds an overview and some example commands of how GeneLab processes Affymetrix microarray datasets. Exact processing commands for specific datasets that have been released is available in this repository [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory and is also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

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

*	`database` – 	the AnnotationDb object
*	`keys` – the keys to select records for from the database. 
*	`-z` – specifies to zip the output data directory
*	`trimmed_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* trimmed_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* trimmed_multiqc_output/trimmed_multiqc_report.html (multiqc output html summary)
* trimmed_multiqc_output/trimmed_multiqc_data.zip (zipped directory containing multiqc output data)

<br>

---

## Assembly-based processing
### 4. Sample assembly
```
megahit -1 sample-1-R1-trimmed.fastq.gz -2 sample-1-R2-trimmed.fastq.gz \
        -o sample-1-assembly -t 10 --min-contig-length 500 > sample-1-assembly.log 2>&1
```

**Parameter Definitions:**  

*	`-1 and -2` – specifies the input forward and reverse reads (if single-end data, then neither `-1` nor `-2` are used, instead single-end reads are passed to `-r`)

*	`-o` – specifies output directory

*	`-t` – specifies the number of threads to use

*	`--min-contig-length` – specifies the minimum contig length to write out

*	`> sample-1-assembly.log 2>&1` – sends stdout/stderr to log file


**Input data:**

* *fastq.gz (filtered/trimmed reads)

**Output data:**

* sample-1-assembly/final.contigs.fa (assembly file)
* sample-1-assembly.log (log file)

<br>

---

### 5. Renaming contigs and summarizing assemblies

#### 5a. Renaming contig headers
```
bit-rename-fasta-headers -i sample-1-assembly/final.contigs.fa -w c_sample-1 -o sample-1-assembly.fasta
```

**Parameter Definitions:**  

*	`-i` – input fasta file

*	`-w` – wanted header prefix (a number will be appended for each contig), starts with a “c_” to ensure they won’t start with a number which can be problematic

*	`-o` – output fasta file


**Input data:**

* sample-1-assembly/final.contigs.fa (assembly file)

**Output files:**

* sample-1-assembly.fasta (contig-renamed assembly file)


#### 5b. Summarizing assemblies

```
bit-summarize-assembly -o assembly-summaries.tsv *assembly.fasta
```

**Parameter Definitions:**  

*	`-o` – output summary table

*	– multiple input assemblies can be provided as positional arguments


**Input data:**

* *-assembly.fasta (contig-renamed assembly files)

**Output files:**

* assembly-summaries.tsv (table of assembly summary statistics)

<br>

---

### 6. Gene prediction
```
prodigal -a sample-1-genes.faa -d sample-1-genes.fasta -f gff -p meta -c -q \
         -o sample-1-genes.gff -i sample-1-assembly.fasta
```
**Parameter Definitions:**

*	`-a` – specifies the output amino acid sequences file

*	`-d` – specifies the output nucleotide sequences file

*	`-f` – specifies the output format gene-calls file

*	`-p` – specifies which mode to run the gene-caller in 

*	`-c` – no incomplete genes reported 

*	`-q` – run in quiet mode (don’t output process on each contig) 

*	`-o` – specifies the name of the output gene-calls file 

*	`-i` – specifies the input assembly

**Input data:**

* sample-1-assembly.fasta (assembly file)

**Output data:**

* sample-1-genes.faa (gene-calls amino-acid fasta file)
* sample-1-genes.fasta (gene-calls nucleotide fasta file)
* sample-1-genes.gff (gene-calls in general feature format)

<br>

---

### 7. Functional annotation
> **Notes**  
> The annotation process overwrites the same temporary directory by default. So if running multiple processses at a time, it is necessary to specify a specific temporary directory with the `--tmp-dir` argument as shown below.


#### 7a. Downloading reference database of HMM models (only needs to be done once)

```
curl -LO ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
curl -LO ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar -xzvf profiles.tar.gz
gunzip ko_list.gz 
```

#### 7b. Running KEGG annotation
```
exec_annotation -p profiles/ -k ko_list --cpu 15 -f detail-tsv -o sample-1-KO-tab.tmp \
                --tmp-dir sample-1-tmp-KO --report-unannotated sample-1-genes.faa 
```

**Parameter Definitions:**
*	`-p` – specifies the directory holding the downloaded reference HMMs

*	`-k` – specifies the downloaded reference KO  (Kegg Orthology) terms 

*	`--cpu` – specifies the number of searches to run in parallel

*	`-f` – specifies the output format

*	`-o` – specifies the output file name

*	`--tmp-dir` – specifies the temporary directory to write to (needed if running more than one process concurrently, see Notes above)

*	`--report-unannotated` – specifies to generate an output for each entry

*	`sample-1-genes.faa` – the input file is specified as a positional argument


**Input data:**

* sample-1-genes.faa (amino-acid fasta file)
* profiles/ (reference directory holding the KO HMMs)
* ko_list (reference list of KOs to scan for)

**Output data:**

* sample-1-KO-tab.tmp (table of KO annotations assigned to gene IDs)


#### 7c. Filtering output to retain only those passing the KO-specific score and top hits
```
bit-filter-KOFamScan-results -i sample-1-KO-tab.tmp -o sample-1-annotations.tsv

  # removing temporary files
rm -rf sample-1-tmp-KO/ sample-1-KO-annots.tmp
```

**Parameter Definitions:**  

*	`-i` – specifies the input table

*	`-o` – specifies the output table


**Input data:**

* sample-1-KO-tab.tmp (table of KO annotations assigned to gene IDs)

**Output data:**

* sample-1-annotations.tsv (table of KO annotations assigned to gene IDs)

<br>

---

### 8. Taxonomic classification

#### 8a. Pulling and un-packing pre-built reference db (only needs to be done once)
```
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz
tar -xvzf CAT_prepare_20200618.tar.gz
```

#### 8b. Running taxonomic classification
```
CAT contigs -c sample-1-assembly.fasta -d CAT_prepare_20200618/2020-06-18_database/ \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ -p sample-1-genes.faa \
            -o sample-1-tax-out.tmp -n 15 -r 3 --top 3 --I_know_what_Im_doing
```

**Parameter Definitions:**  

*	`-c` – specifies the input assembly fasta file

*	`-d` – specifies the CAT reference sequence database

*	`-t` – specifies the CAT reference taxonomy database

*	`-p` – specifies the input protein fasta file

*	`-o` – specifies the output prefix

*	`-n` – specifies the number of cores to use

*	`-r` – specifies the number of top protein hits to consider in assigning tax

*	`--top` – specifies the number of protein alignments to store

*	`--I_know_what_Im_doing` – allows us to alter the `--top` parameter


**Input data:**

* sample-1-assembly.fasta (assembly file)
* sample-1-genes.faa (gene-calls amino-acid fasta file)

**Output data:**

* sample-1-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)
* sample-1-tax-out.tmp.contig2classification.txt (contig taxonomy file)

#### 8c. Adding taxonomy info from taxids to genes
```
CAT add_names -i sample-1-tax-out.tmp.ORF2LCA.txt -o sample-1-gene-tax-out.tmp \
              -t CAT_prepare_20200618/2020-06-18_taxonomy/ --only_official
```

**Parameter Definitions:**  

*	`-i` – specifies the input taxonomy file

*	`-o` – specifies the output file 

*	`-t` – specifies the CAT reference taxonomy database

*	`--only_official` – specifies to add only standard taxonomic ranks

**Input data:**

* sample-1-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)

**Output data:**

* sample-1-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)



#### 8d. Adding taxonomy info from taxids to contigs
```
CAT add_names -i sample-1-tax-out.tmp.contig2classification.txt -o sample-1-contig-tax-out.tmp \
              -t CAT-ref/2020-06-18_taxonomy/ --only_official
```

**Parameter Definitions:**  

*	`-i` – specifies the input taxonomy file

*	`-o` – specifies the output file 

*	`-t` – specifies the CAT reference taxonomy database

*	`--only_official` – specifies to add only standard taxonomic ranks


**Input data:**

* sample-1-tax-out.tmp.contig2classification.txt (contig taxonomy file)

**Output data:**

* sample-1-contig-tax-out.tmp (contig taxonomy file with lineage info added)


#### 8e. Formatting gene-level output with awk and sed
```
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $2 == "lineage" ) { print $1,$2,$4,$5,$6,$7,$8,$9,$10 } \
    else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) \
    { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } else { n=split($2,lineage,";"); \
    print $1,lineage[n],$4,$5,$6,$7,$8,$9,$10 } } ' sample-1-gene-tax-out.tmp | \
    sed 's/not classified/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# ORF/gene_ID/' | \
    sed 's/lineage/taxid/' | sed 's/\*//g' > sample-1-gene-tax-out.tsv
```

#### 8f. Formatting contig-level output with awk and sed
```
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $2 == "classification" ) { print $1,$4,$6,$7,$8,$9,$10,$11,$12 } \
    else if ( $2 == "unclassified" ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
    else { n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 } } ' sample-1-contig-tax-out.tmp | \
    sed 's/not classified/NA/g' | sed 's/superkingdom/domain/' | sed 's/: [0-9\.]*//g' | sed 's/^# contig/contig_ID/' | \
    sed 's/lineage/taxid/' | sed 's/\*//g' > sample-1-contig-tax-out.tsv

  # clearing intermediate files
rm sample-1*.tmp*
```

**Input data:**

* sample-1-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)
* sample-1-contig-tax-out.tmp (contig taxonomy file with lineage info added)


**Output data:**

* sample-1-gene-tax-out.tsv (gene-calls taxonomy file with lineage info added reformatted)
* sample-1-contig-tax-out.tsv (contig taxonomy file with lineage info added reformatted)

<br>

---

### 9. Read-mapping

#### 9a. Building reference index
```
bowtie2-build sample-1-assembly.fasta sample-1-assembly-bt-index
```

**Parameter Definitions:**  

*	`sample-1-assembly.fasta` - first positional argument specifies the input assembly

*	`sample-1-assembly-bt-index` - second positional argument specifies the prefix of the output index files


#### 9b. Performing mapping, conversion to bam, and sorting
```
bowtie2 --threads 15 -x sample-1-assembly-bt-index -1 sample-1-R1-trimmed.fastq.gz \
        -2 sample-1-R2-trimmed.fastq.gz 2> sample-1-mapping.log | samtools view -b | samtools sort -@ 15 > sample-1.bam
```

**Parameter Definitions:**  

*	`--threads` – specifies the number of threads to run in parallel

*	`-x` – specifies the prefix of the reference index files to map to (generated in the previous `bowtie2-build` step

*	`-1 and -2` – specifies the forward and reverse reads to map (if single-end data, neither `-1` nor `-2` are provided, and the single-end reads are passed to `-r`)

* `2> sample-1-mapping.log` – capture the printed summary results in a log file

*	`samtools view -b` – convert the output directly to bam format (compressed)

*	`samtools sort -@` – sort the bam file using the specified number of threads

*	`>` – redirect the output to a file

#### 9c. Indexing
```
samtools index -@ 15 sample-1.bam 
```

**Parameter Definitions:**  
*	`-@` – set number of threads to use 

*	`sample-1.bam` - input bam file is provided as a positional argument as generated from the above mapping step

**Input data:**

* sample-1-assembly.fasta (assembly file)
* *.fastq.gz (filtered/trimmed reads)

**Output data:**

* sample-1.bam (mapping file)
* sample-1.bam.bai (bam index file)
* sample-1-mapping.log (read-mapping log file)

<br>

---

### 10. Getting coverage information and filtering based on detection
> **Notes**  
> “Detection” is a metric of what proportion of a reference sequence recruited reads (see [here](http://merenlab.org/2017/05/08/anvio-views/#detection)). Filtering based on detection is one way of helping to mitigate non-specific read-recruitment.

#### 10a. Filtering coverage levels based on detection

```
  # pileup.sh comes from the bbduk.sh package
pileup.sh -in sample-1.bam fastaorf=sample-1-genes.fasta outorf=sample-1-gene-cov-and-det.tmp \
          out=sample-1-contig-cov-and-det.tmp
```

**Parameter Definitions:**  

*	`-in` – the input bam file

*	`fastaorf=` – input gene-calls nucleotide fasta file

*	`outorf=` – the output gene-coverage tsv file

*	`out=` – the output contig-coverage tsv file


#### 10b. Filtering gene coverage based on requiring 50% detection and parsing down to just gene ID and coverage
```
grep -v "#" sample-1-gene-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $10 <= 0.5 ) $4 = 0 } \
     { print $1,$4 } ' > sample-1-gene-cov.tmp

cat <( printf "gene_ID\tcoverage\n" ) sample-1-gene-cov.tmp > sample-1-gene-coverages.tsv
```

Filtering contig coverage based on requiring 50% detection and parsing down to just contig ID and coverage:
```
grep -v "#" sample-1-contig-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $5 <= 50 ) $2 = 0 } \
     { print $1,$2 } ' > sample-1-contig-cov.tmp

cat <( printf "contig_ID\tcoverage\n" ) sample-1-contig-cov.tmp > sample-1-contig-coverages.tsv

  # removing intermediate files

rm sample-1-*.tmp
```

**Input data:**

* sample-1.bam (mapping file)
* sample-1-genes.fasta (gene-calls nucleotide fasta file)

**Output data:**

* sample-1-gene-coverages.tsv (table with gene-level coverages)
* sample-1-contig-coverages.tsv (table with contig-level coverages)

<br>

---

### 11. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk`, all are standard in any Unix-like environment.  

```
paste <( tail -n +2 sample-1-gene-coverages.tsv | sort -V -k 1 ) <( tail -n +2 sample-1-annotations.tsv | sort -V -k 1 | cut -f 2- ) \
      <( tail -n +2 sample-1-gene-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-1-gene-tab.tmp

paste <( head -n 1 sample-1-gene-coverages.tsv ) <( head -n 1 sample-1-annotations.tsv | cut -f 2- ) \
      <( head -n 1 sample-1-gene-tax-out.tsv | cut -f 2- ) > sample-1-header.tmp

cat sample-1-header.tmp sample-1-gene-tab.tmp > sample-1-gene-coverage-annotation-and-tax.tsv

  # removing intermediate files
rm sample-1*tmp sample-1-gene-coverages.tsv sample-1-annotations.tsv sample-1-gene-tax-out.tsv
```

**Input data:**

* sample-1-gene-coverages.tsv (table with gene-level coverages from step 10)
* sample-1-annotations.tsv (table of KO annotations assigned to gene IDs from step 7)
* sample-1-gene-tax-out.tsv (gene-level taxonomic classifications from step 8)


**Output data:**

* sample-1-gene-coverage-annotation-and-tax.tsv (table with combined gene coverage, annotation, and taxonomy info)

<br>

---

### 12. Combining contig-level coverage and taxonomy into one table for each sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk`, all are standard in any Unix-like environment.  

```
paste <( tail -n +2 sample-1-contig-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-1-contig-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-1-contig.tmp

paste <( head -n 1 sample-1-contig-coverages.tsv ) <( head -n 1 sample-1-contig-tax-out.tsv | cut -f 2- ) \
      > sample-1-contig-header.tmp
      
cat sample-1-contig-header.tmp sample-1-contig.tmp > sample-1-contig-coverage-and-tax.tsv

  # removing intermediate files
rm sample-1*tmp sample-1-contig-coverages.tsv sample-1-contig-tax-out.tsv
```

**Input data:**

* sample-1-contig-coverages.tsv (table with contig-level coverages from step 10)
* sample-1-contig-tax-out.tsv (contig-level taxonomic classifications from step 8)


**Output data:**

* sample-1-contig-coverage-and-tax.tsv (table with combined contig coverage and taxonomy info)

<br>

---

### 13. Generating normalized, gene-level-coverage summary tables of KO-annotations and taxonomy across samples
> **Notes**  
> * To combine across samples to generate these summary tables, we need the same "units". This is done for annotations based on the assigned KO terms, and all non-annotated functions are included together as "Not annotated". It is done for taxonomic classifications based on taxids (full lineages included in the table), and any not classified are included together as "Not classified". 
> * The values we are working with are coverage per gene (so they are number of bps recruited to the gene normalized by the length of the gene). These have been normalized by making the total coverage of a sample 1,000,000 and setting each individual gene-level coverage its proportion of that 1,000,000 total. So basically percent, but out of 1,000,000 instead of 100 to make the numbers more friendly. 

```
bit-GL-combine-KO-and-tax-tables *-gene-coverage-annotation-and-tax.tsv -o GLDS-286
```

**Parameter Definitions:**  

*	takes positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above

-	`-o` – specifies the output prefix (e.g. as above, will generate “GLDS-286-KO-function-coverages.tsv” and “GLDS-286-taxonomy-coverages.tsv”


**Input data:**

* *-gene-coverage-annotation-and-tax.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples from step 12)

**Output data:**

* GLDS-286-KO-function-coverages.tsv (table with all samples combined based on KO annotations; normalized to coverage per million genes covered)
* GLDS-286-taxonomy-coverages.tsv (table with all samples combined based on gene-level taxonomic classifications; normalized to coverage per million genes covered)

<br>

---

### 14. **M**etagenome-**A**ssembled **G**enome (MAG) recovery

#### 14a. Binning contigs
```
jgi_summarize_bam_contig_depths --outputDepth sample-1-depth.tsv --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta sample-1-assembly.fasta sample-1.bam

metabat2  --inFile sample-1-assembly.fasta --outFile sample-1 --abdFile sample-1-depth.tsv -t 4
```

**Parameter Definitions:**  

*  `--outputDepth` – specifies the output depth file
*  `--percentIdentity` – minimum end-to-end percent identity of a mapped read to be included
*  `--minContigLength` – minimum contig length to include
*  `--minContigDepth` – minimum contig depth to include
*  `--referenceFasta` – the assembly fasta file generated in step 4
*  `sample-1.bam` – final positional arguments are the bam files generated in step 9
*  `--inFile` - the assembly fasta file generated in step 4
*  `--outFile` - the prefix of the identified bins output files
*  `--abdFile` - the depth file generated by the previous `jgi_summarize_bam_contig_depths` command
*  `-t` - specifies number of threads to use


**Input data:**

* sample-1-assembly.fasta (assembly fasta file created in step 4)
* sample-1.bam (bam file created in step 9)

**Output data:**

* sample-1-depth.tsv (tab-delimited summary of coverages)
* sample-1-bin\*.fa (fasta files of recovered bins)

#### 14b. Bin quality assessment
Utilizes the default `checkm` database available [here](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz), `checkm_data_2015_01_16.tar.gz`.

```
checkm lineage_wf -f checkm-bins-summary.tsv --tab_table -x fa ./ checkm-output-dir
```

**Parameter Definitions:**  

*  `lineage_wf` – specifies the workflow being utilized
*  `-f` – specifies the output summary file
*  `--tab_table` – specifies the output summary file should be a tab-delimited table
*  `-x` – specifies the extension that is on the bin fasta files that are being assessed
*  `./` – first positional argument at end specifies the directory holding the bins generated in step 14a
*  `checkm-output-dir` – second positional argument at end specifies the primary checkm output directory with detailed information

**Input data:**

* bin fasta files generated by step 14a

**Output data:**

* checkm-bins-summary.tsv (tab-delimited file with quality estimates per bin)
* checkm-output-dir (directory holding detailed checkm outputs)

#### 14c. Filtering MAGs

```
cat <( head -n 1 checkm-bins-summary.tsv ) \
    <( awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' checkm-bins-summary.tsv | sed 's/bin./MAG-/' ) \
    > checkm-MAGs-summary.tsv
    
# copying bins into a MAGs directory in order to run tax classification
awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' checkm-bins-summary.tsv | cut -f 1 > MAG-bin-IDs.tmp

mkdir MAGs
for ID in MAG-bin-IDs.tmp
do
    MAG_ID=$(echo $ID | sed 's/bin./MAG-/')
    cp ${ID}.fa MAGs/${MAG_ID}.fa
done
```

**Input data:**

* checkm-bins-summary.tsv (tab-delimited file with quality estimates per bin)

**Output data:**

* checkm-MAGs-summary.tsv (tab-delimited file with quality estimates per MAG)
* MAGs/\*.fa (directory holding high-quality MAGs)



#### 14d. MAG taxonomic classification
Uses default `gtdbtk` database setup with program's `download.sh` command.

```
gtdbtk classify_wf --genome_dir MAGs/ -x fa --out_dir gtdbtk-output-dir 
```

**Parameter Definitions:**  

*  `classify_wf` – specifies the workflow being utilized
*  `--genome_dir` – specifies the directory holding the MAGs generated in step 14c
*  `-x` – specifies the extension that is on the MAG fasta files that are being taxonomically classified
*  `-out_dir` – specifies the output directory

**Input data:**

* MAGs/\*.fa (directory holding high-quality MAGs)

**Output data:**

* gtdbtk-output-dir/gtdbtk.\*.summary.tsv (files with assigned taxonomy and info)

<br>

---

## Read-based processing
### 15. Taxonomic and functional profiling
The following uses the `humann3` and `metaphlan3` reference databases downloaded on 26-Sept-2020 as follows:

```bash
humann_databases --download chocophlan full
humann_databases --download uniref uniref90_diamond 
humann_databases --download utility_mapping full 
metaphlan --install
```

#### 15a. Running humann3 (which also runs metaphlan3)
```bash
  # forward and reverse reads need to be provided combined if paired-end (if not paired-end, single-end reads are provided to the --input argument next)
cat sample-1-R1-trimmed.fastq.gz sample-1-R2-trimmed.fastq.gz > sample-1-combined.fastq.gz

humann --input sample-1-combined.fastq.gz --output sample-1-humann3-out-dir --threads 15 \
       --output-basename sample-1 --metaphlan-options "--unknown_estimation --add_viruses \
       --sample_id sample-1"
```

**Parameter Definitions:**  

*	`--input` – specifies the input combined forward and reverse reads (if paired-end)

*	`--output` – specifies output directory

*	`--threads` – specifies the number of threads to use

*	`--output-basename` – specifies prefix of the output files

*	`--metaphlan-options` – options to be passed to metaphlan
	* `--unknown_estimation` – include unclassified in estimated relative abundances
	* `--add_viruses` – include viruses in the reference database
	* `--sample_id` – specifies the sample identifier we want in the table (rather than full filename)


#### 15b. Merging multiple sample functional profiles into one table
```bash
  # they need to be in their own directories
mkdir genefamily-results/ pathabundance-results/ pathcoverage-results/

  # copying results from previous running humann3 step (14a) to get them all together in their own directories
cp *-humann3-out-dir/*genefamilies.tsv genefamily-results/
cp *-humann3-out-dir/*abundance.tsv pathabundance-results/
cp *-humann3-out-dir/*coverage.tsv pathcoverage-results/

humann_join_tables -i genefamily-results/ -o gene-families.tsv
humann_join_tables -i pathabundance-results/ -o path-abundances.tsv
humann_join_tables -i pathcoverage-results/ -o path-coverages.tsv
```

**Parameter Definitions:**  

*	`-i` – the directory holding the input tables

*	`-o` – the name of the output combined table


#### 15c. Splitting results tables
The read-based functional annotation tables have taxonomic info and non-taxonomic info mixed together initially. `humann` comes with a helper script to split these. Here we are using that to generate both non-taxonomically grouped functional info files and taxonomically grouped ones.

```bash
humann_split_stratified_table -i gene-families.tsv -o ./
mv gene-families_stratified.tsv gene-families-grouped-by-taxa.tsv
mv gene-families_unstratified.tsv gene-families.tsv

humann_split_stratified_table -i path-abundances.tsv -o ./
mv path-abundances_stratified.tsv path-abundances-grouped-by-taxa.tsv
mv path-abundances_unstratified.tsv path-abundances.tsv

humann2_split_stratified_table -i path-coverages.tsv -o ./
mv path-coverages_stratified.tsv path-coverages-grouped-by-taxa.tsv
mv path-coverages_unstratified.tsv path-coverages.tsv
```

**Parameter Definitions:**  

*	`-i` – the input combined table

*	`-o` – output directory (here specifying current directory)


#### 15d. Normalizing gene families and pathway abundance tables
This generates some normalized tables of the read-based functional outputs from humann that are more readily suitable for across sample comparisons.

```bash
humann_renorm_table -i gene-families.tsv -o gene-families-cpm.tsv --update-snames
humann_renorm_table -i path-abundances.tsv -o path-abundances-cpm.tsv --update-snames
```

**Parameter Definitions:**  

*	`-i` – the input combined table

*	`-o` – name of the output normalized table

*	`--update-snames` – change suffix of column names in tables to "-CPM"


#### 15e. Generating a normalized gene-family table that is grouped by Kegg Orthologs (KOs)

```bash
humann_regroup_table -i gene-families.tsv -g uniref90_ko | humann_rename_table -n kegg-orthology | \
                     humann_renorm_table -o gene-families-KO-cpm.tsv --update-snames
```

**Parameter Definitions:**  

*	`-i` – the input table

*	`-g` – the map to use to group uniref IDs into Kegg Orthologs

*	`|` – sending that output into the next humann command to add human-readable Kegg Orthology names

*	`-n` – specifying we are converting Kegg orthology IDs into Kegg orthology human-readable names

*	`|` – sending that output into the next humann command to normalize to copies-per-million

*	`-o` – specifying the final output file name

*  `--update-snames` – change suffix of column names in tables to "-CPM"

#### 15f. Combining taxonomy tables

```bash
merge_metaphlan_tables.py *-humann3-out-dir/*_humann_temp/*_metaphlan_bugs_list.tsv > metaphlan-taxonomy.tsv
```

**Parameter Definitions:**  

*	input metaphlan tables are provided as position arguments (produced during humann3 run above, step 14a)

*  `>` – output is redirected from stdout to a file


**Input data:**

* *fastq.gz (filtered/trimmed reads from step 2, forward and reverse reads concatenated if paired-end)

**Output data:**

* gene-families.tsv (gene-family abundances) 
* gene-families-grouped-by-taxa.tsv (gene-family abundances grouped by taxa)
* gene-families-cpm.tsv (gene-family abundances normalized to copies-per-million)
* gene-families-KO-cpm.tsv (KO term abundances normalized to copies-per-million)
* pathway-abundances.tsv (pathway abundances)
* pathway-abundances-grouped-by-taxa.tsv (pathway abundances grouped by taxa)
* pathway-abundances-cpm.tsv (pathway abundances normalized to copies-per-million)
* pathway-coverages.tsv (pathway coverages)
* pathway-coverages-grouped-by-taxa.tsv (pathway coverages grouped by taxa)
* metaphlan-taxonomy.tsv (metaphlan estimated taxonomic relative abundances)

---
