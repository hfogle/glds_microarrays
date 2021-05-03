#!/usr/bin/env Rscript


library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Dataset file paths", metavar="character"),
  make_option(c("-i", "--isa"), type="character", default=NULL, 
              help="Study ISAtab.zip file path", metavar="character"),
  make_option(c("-p", "--probe"), type="character", default=NULL, 
              help="Probe annotation file path", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="Specify organism species", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="Organism annotation file path", metavar="character"),
  make_option(c("-g", "--glds"), type="character", default=NULL, 
              help="GeneLab study identifier", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Local input directory", metavar="character"),
  make_option(c("-m", "--platform"), type="character", default=NULL, 
              help="Microarray Platform", metavar="character"), 
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output directory path [default= %default]", metavar="character"),
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

### Determine Platform

### Selects appropriate package for plots
if ((opt$platform == "Affymetrix Expression") | (opt$platform == "Affymetrix ST") | (opt$platform == "NimbleGen")){
  plots <- "OLIGO"
}else if((opt$platform == "Illumina BeadChip") | (opt$platform == "Agilent") | (opt$platform == "Agilent GenePix")){
  plots <- "LIMMA"
}else {
  plots <- NULL
}

### Parse ISA

source("microarray_functions.R")
all_targets<-buildTargets(opt)
if (length(all_targets) == 1){
  targets <- all_targets[[1]]
}else{
  targets <- all_targets[[which_assay]]
}

### Select Filters (MaxIQR, Annotated, Control Probes, Low Expression)

opt$model <- "Group Contrast"
opt$estimation <- "Max IQR"

opt$filter_nonexpressed <- TRUE
opt$filter_control <- TRUE
opt$filter_nonannotated <- TRUE

### Call appropriate platform processing script

if(opt$platform=="affymetrix"){
  source(file.path("platform","affymetrix.R"))
}else if(opt$platform=="agilent_one_channel"){
  source(file.path("platform","agilent_one_channel.R"))
}else if(opt$platform=="agilent_two_channel"){
  source(file.path("platform","agilent_two_channel.R"))
}else{
  print_help(opt_parser)
  stop("Platform not currently supported.n", call.=FALSE)
}
