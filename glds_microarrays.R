#!/usr/bin/env Rscript


library("optparse")

option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="Dataset file paths", metavar="character"),
  make_option(c("-i", "--isa"), type="character", default=NULL, 
              help="Study ISAtab.zip file path", metavar="character"),
  make_option(c("-p", "--probe"), type="character", default=NULL, 
              help="Probe annotation file path", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="Specify organism species (e.g. Arabidopsis thaliana)", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="Organism annotation file path", metavar="character"),
  make_option(c("-g", "--glds"), type="character", default=NULL, 
              help="GeneLab study identifier (e.g. GLDS-121)", metavar="character"),
  make_option(c("-m", "--platform"), type="character", default=NULL, 
              help="Microarray Platform", metavar="character"), 
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Curated input directory path", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=getwd(), 
              help="output parent directory path [default= %default]", metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$glds)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (GLDS)", call.=FALSE)
}
### Get organism annotation package
options(connectionObserver = NULL)
organism_table <- read.csv(file = file.path(getwd(),"organisms.csv"), header = TRUE, stringsAsFactors = FALSE)
ann.dbi <- organism_table$annotations[organism_table$species == opt$species] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
  BiocManager::install(ann.dbi, ask = FALSE)
  library(ann.dbi, character.only=TRUE)
}
#require(ann.dbi, character.only = TRUE, quietly = TRUE)
cat("\nOrganism annotation set loaded: ",ann.dbi,"\n")


### Get Study Files
if (is.null(opt$files)){
  cat("\nStart processing ",opt$glds, "\n")
}

if (!is.null(opt$dir)){
  cat("Building from local curated directory\n")
  opt$isa <- list.files(file.path(getwd(),opt$dir,"Metadata"),pattern = "*ISA.zip", full.names = TRUE)
  opt$probe <- list.files(file.path(getwd(),opt$dir,"Metadata"),pattern = "*annotation*", full.names = TRUE)
  opt$files <- list.files(file.path(getwd(),opt$dir,"00-RawData"),full.names = TRUE)
  str(opt)
  
  tempin <- tempdir()
  unlink(list.files(tempin, full.names = TRUE))
  dir.create(file.path(tempin,"00-RawData"), showWarnings = FALSE)
  file.copy(from = opt$files, to = file.path(tempin,"00-RawData"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
  opt$files <- list.files(file.path(tempin,"00-RawData"),full.names = TRUE)
  dir.create(file.path(tempin,"Metadata"),showWarnings = FALSE)
  file.copy(from = opt$isa, to = file.path(tempin,"Metadata"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
  opt$isa <- list.files(file.path(tempin,"Metadata"),pattern = "*ISA.zip", full.names = TRUE)
  file.copy(from = opt$probe, to = file.path(tempin,"Metadata"), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
  opt$probe <- list.files(file.path(tempin,"Metadata"),pattern = "*annotation*", full.names = TRUE)
}


### Determine Platform

### Selects appropriate package for plots
plots <- NULL
if ((opt$platform == "Affymetrix") | (opt$platform == "NimbleGen")){
  plots <- "OLIGO"
}else if((opt$platform == "Illumina BeadChip") | (opt$platform == "Agilent") | (opt$platform == "Agilent GenePix")){
  plots <- "LIMMA"
}else {
  plots <- NULL
}
cat("\nPlots generate by package ",plots,"\n")

### Parse ISA

source("microarray_functions.R")
all_targets<-buildTargets(opt)
if (length(all_targets) == 1){
  targets <- all_targets[[1]]
}else{
  targets <- all_targets[[which_assay]]
}
cat("\nTargets\n")
str(targets)



### Call appropriate platform processing script

if(opt$platform=="Affymetrix"){
  source(file.path(getwd(),"platforms","affymetrix.R"))
}else if(opt$platform=="agilent_one_channel"){
  source(file.path(getwd(),"platforms","agilent_one_channel.R"))
}else if(opt$platform=="agilent_two_channel"){
  source(file.path(getwd(),"platforms","agilent_two_channel.R"))
}else{
  print_help(opt_parser)
  stop("Platform not currently supported.n", call.=FALSE)
}
