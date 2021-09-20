# GeneLab Two-Channel Expression Microarray Data Processing Pipeline

This is an R Markdown Interactive Document for processing GeneLab curated datasets that contain two-color gene expression assays. It can be run within the RStudio IDE with the Knit HTML button or on an R console with the command:  
rmarkdown::render("glds_two_channel_arrays.Rmd")

An html report file and dynamic site will be generated from parameters defined in the manual_entry code block. Interactive parameter selection is available using the Knit with Parameters option. Processed data files will also be exported to the working directory.  

Details about design model and normalization function parameters can be found in the Limma package documentation.[link](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)  

***
### Input Parameters

* **GLDS Accession #** --- GeneLab repository accession number of the dataset or any string to create a seperate output folder (e.g. GLDS-8, GLDS8, Test1)
* **Organism** --- Supported species annotations are listed in the YAML header
* **Microarray Data Format** --- The Platform or Scanner file formatting of the raw data. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'image_aquisition' 'TYPE'
  + Agilent GenePix --- Agilent microarray platform imaged with GenePix software
  + Agilent --- Agilent microarray platform imaged with Agilent software
  + GenePix --- GenePix software formatting and intensity estimation
  + Spot --- SPOT software formatting with mean foreground estimation and morph background estimation
  + Imagene --- Imagene software formatting and intensity estimation
  + Bluefuse --- Bluefuse software formatting and intensity estimation
  + Generic --- Generic spot array
* **Two Channel Configuration** --- The sample to color channel configuration of the dataset. Dye-swapping arrangements are detected.
  + Replicate Array --- Only two samples present across multiple arrays
  + Common Reference --- A pool reference sample is present on each array
  + Direct Two-Color --- Samples are compared directly by competitive hybridization on the same arrays
  + Separate Channels ---  Converts a two-color experiment into a single channel experiment with twice as many arrays but with a technical pairing between the two channels
* **Background Correction** --- Options supported by limma. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'feature_extraction' 'TYPE'
  + normexp --- is recommended as the default background correction method.
  + subtract ---
  + morph ---
* **Within Array Normalization** --- Options supported by limma. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'feature_extraction' 'TYPE'
  + loess --- Loess within array normalization recommended for Agilent arrays
  + printtiploess --- Print-tip recommended in general for others
  + robustspline ---
* **Between Array Normalization** --- Options supported by limma. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'feature_extraction' 'TYPE'
  + Aquantile --- recommended as default. Applies quantile normalization to A-values between arrays
  + quantile --- normalizes directly to the individual red and green intensities
* **Primary Annotation Keytype** --- Gene identifier type for probe annotation from [Ann.dbi](https://www.bioconductor.org/packages/release/data/annotation/) sources
* **Filters** --- Probe level filtering options to indicate which probes to remove prior to differential gene expression (DGE) analysis. By default, Non-Expressed, Control, and Non-Annotated Probes are all removed before DGE analysis to improve statistics and avoid modeling errors.
* **Raw Data Files** --- Select the raw data files to be analyzed
* **Probe Annotation File** --- Select the Annotation file provided by the dataset submitter. This file can be downloaded from the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for each microarray dataset under 'STUDY FILES' -> 'Microarray Data Files' -> &ast;adf.txt (from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)) or GPL*.soft (from [GEO](https://www.ncbi.nlm.nih.gov/geo/)), or a platform specific file supplied by the dataset submitter.
* **ISA Metadata File** --- GeneLab ISAtab study metadata file, which can be downloaded from the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files' -> *ISA.zip

***
### Required R Packages available from the CRAN repository:  

Package   | Documentation
----------|--------------
shiny     | [https://cran.r-project.org/web/packages/shiny/index.html](https://cran.r-project.org/web/packages/shiny/index.html)
rmarkdown | [https://cran.r-project.org/web/packages/rmarkdown/index.html](https://cran.r-project.org/web/packages/rmarkdown/index.html)
knitr     | [https://cran.r-project.org/web/packages/knitr/index.html](https://cran.r-project.org/web/packages/knitr/index.html)
xfun      | [https://cran.r-project.org/web/packages/xfun/index.html](https://cran.r-project.org/web/packages/xfun/index.html)
DT        | [https://cran.r-project.org/web/packages/DT/index.html](https://cran.r-project.org/web/packages/DT/index.html)
R.utils   | [https://cran.r-project.org/web/packages/R.utils/index.html](https://cran.r-project.org/web/packages/R.utils/index.html)
dplyr     | [https://cran.r-project.org/web/packages/dplyr/index.html](https://cran.r-project.org/web/packages/dplyr/index.html)
tydyr     | [https://cran.r-project.org/web/packages/tidyr/index.html](https://cran.r-project.org/web/packages/tidyr/index.html)
statmod   | [https://cran.r-project.org/web/packages/statmod/index.html](https://cran.r-project.org/web/packages/statmod/index.html)


### Required R Packages available from the Bioconductor repository:

Package       | Documentation
--------------|--------------
Risa          | [https://www.bioconductor.org/packages/release/bioc/html/Risa.html](https://www.bioconductor.org/packages/release/bioc/html/Risa.html)
limma         | [https://bioconductor.org/packages/release/bioc/html/limma.html](https://bioconductor.org/packages/release/bioc/html/limma.html)
GEOquery      | [https://bioconductor.org/packages/release/bioc/html/GEOquery.html](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)
ArrayExpress  | [https://www.bioconductor.org/packages/release/bioc/html/ArrayExpress.html](https://www.bioconductor.org/packages/release/bioc/html/ArrayExpress.html)
STRINGdb      | [https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)
AnnotationDbi | [https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html](https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)

