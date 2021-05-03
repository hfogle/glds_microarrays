### Agilent Two Channel Microarray Processing

### Import Raw Data

library(limma)
path <- dirname(opt$datafiles$datapath[1])
files_in_assay <- sapply(targets$t1$FileName,function(x) grep(x,opt$datafiles$name))
files_in_assay <- opt$datafiles$name[files_in_assay]

if (opt$format=="GPR"){
  RG <- read.maimages(files_in_assay, source="genepix", path=path,wt.fun=wtflags(weight=0,cutoff=-50))
}else if (opt$format=="TXT"){
  RG <- read.maimages(files_in_assay, source="agilent", path=path)
}

dir.create(file.path(workdir,"Processed_Data"), showWarnings = FALSE)
dir.create(file.path(workdir,"Processed_Data",opt$glds), showWarnings = FALSE)
dir.create(file.path(workdir,"Processed_Data",opt$glds,"00-RawData"), showWarnings = FALSE)
path <- file.path(workdir,"Processed_Data",opt$glds,"00-RawData")
file.copy(from = opt$datafiles$datapath, to = file.path(path,opt$datafiles$name), overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
rm(path)


### Import Probe Annotation

if (grepl("\\.soft$", opt$probefile$name) || grepl("\\.gpl.txt$", opt$probefile$name)){
    probedata<- GEOquery:::parseGPL(opt$probefile$datapath)
    probedata<- probedata@dataTable@table
    probedata <- probedata[with(probedata, order(ROW,COL)),] # sort by row & col to match RG order
    RG$genes <- dplyr::left_join(RG$genes,probedata,by = c("Name" = "NAME"))
} 

if (grepl("\\.adf.txt$", opt$probefile$name)){
  probedata<- ArrayExpress:::readFeatures(opt$probefile$name, dirname(opt$probefile$datapath),ArrayExpress:::skipADFheader(opt$probefile$name, dirname(opt$probefile$datapath),proc=F))
  
}

### Import Organism Annotation
keytype = opt$primary
organism_table <- read.csv(file = file.path(getwd(),"organisms.csv"), header = TRUE, stringsAsFactors = FALSE)
ann.dbi <- organism_table$annotations[organism_table$species == opt$organism] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
  BiocManager::install(ann.dbi, ask = FALSE)
  library(ann.dbi, character.only=TRUE)
}

### Map Annotations from Primary Keytype
try(annotation<-data.frame(ID=RG$genes$ID)) #assumes Array ID is primary annotation
try(annotation<-data.frame(ID=RG$genes$GeneName))
try(annotation<-data.frame(ID=RG$genes[,keytype]))
try(annotation<-data.frame(ID=RG$genes$Name)) #assumes Array ID is primary annotation
annotation$ID <- sub("\\.\\d+", "", annotation$ID) #remove any version suffixes from IDs

try(annotation$REFSEQ<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "REFSEQ",multiVals = "first")))
try(annotation$ENSEMBL<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "ENSEMBL",multiVals = "first")))
try(annotation$SYMBOL<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "SYMBOL",multiVals = "first"))) #GENENAME
try(annotation$DESCRIPTION<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "GENENAME",multiVals = "first"))) #DESCRIPTION
try(annotation$TAIR<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "TAIR",multiVals = "first")))
try(annotation$ORF<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "ORF",multiVals = "first")))
try(annotation$ENTREZID<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "ENTREZID",multiVals = "first")))
try(annotation$GOSLIM_ID<-as.character(mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = as.character(annotation$ID),keytype = keytype, column = "GO",multiVals = "first")))

annotation <- data.frame(annotation,stringsAsFactors = FALSE)
is.na(annotation) <- annotation == "NULL"


### Map STRING DB annotations
try({
  string_db <- STRINGdb::STRINGdb$new( version="11", species=organism_table$taxon[organism_table$species == opt$organism],score_threshold=0)
  string_map<-string_db$map(annotation,"REFSEQ",removeUnmappedRows = FALSE, takeFirst = TRUE)
  annotation$STRING_ID <- string_map$STRING_id
  rm(string_map,string_db)
})

annotation.subset <- annotation
RG.annotated <- RG


### Background Correction and Normalization

RG.annotated <- backgroundCorrect(RG.annotated, method="normexp", offset=50,normexp.method="saddle")
MA <- normalizeWithinArrays(RG.annotated, method="loess", weights = RG.annotated$weights)
MA <- normalizeBetweenArrays(MA, method="Aquantile")
MA.summarized <- MA
MA.summarized$genes<-annotation.subset


### Filter out Control Probes
if (opt$filter_control){
  #try(MA.summarized <- MA.summarized[annotation1$CONTROL_TYPE == FALSE,])
  try({
    filter <- which(!(MA.summarized$genes$CONTROLTYPE %in% c("pos","neg")))
    MA.summarized<- MA.summarized[filter,]
    annotation.subset <- annotation.subset[filter,]})
}

### Filter out Non-annotated Probes
if (opt$filter_nonannotated){
  try({
    MA.summarized <- MA.summarized[!is.na(MA.summarized$genes$ENTREZID),]
    annotation.subset <- MA.summarized$genes
    annotation_stats$annotated_genes <- dim(MA.summarized)[1]
  })
}


### Generate Gene Summary Values

MA.summarized <- avereps(MA.summarized, ID = MA.summarized$genes$ID)
annotation.subset <- MA.summarized$genes
annotation_stats$unique_probes <- dim(MA.summarized)[1]

###  Generate Expression Value  Files
dir.create(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"01-NormalizedData"))

#annotation.subset <- dplyr::select(MA.summarized$genes, contains(c("ID")))
colnames(annotation.subset)[1]<-opt$primary
#annotation.subset <-apply(annotation.subset, 1, unlist)
expression <- as.data.frame(MA.summarized$M)
write.table(expression,"normalized.txt",quote=FALSE, append=FALSE, sep = "\t", col.names=NA)

### Write out annotated expression values

expression <- cbind(annotation.subset,expression)
write.table(expression,"normalized-annotated.txt",quote=FALSE, append=FALSE, sep = "\t", row.names = FALSE)

### Write out annotated MA file
save(MA.summarized,file = "normalized-annotated.rda")

### Build Differential Gene Expression Model
setwd(file.path(workdir,"Processed_Data",opt$glds,"00-RawData","QC_Repports"))
library(statmod)

if (opt$model == "Replicate Array"){
 
    uu<-uniqueTargets(targets$t2)
    uuu<-uniqueTargets(targets$t3)
    fff <- factor(targets$t3$Cy3, levels=uuu)
    fit <- lmFit(MA.summarized, ref=uuu[1])
    contrasts<-c(paste(uuu[1],uuu[2],sep = "-"),paste(uuu[2],uuu[1],sep = "-"))
    fit <- eBayes(fit)
    contrast.names<-c(paste(uu[1],uu[2],sep = "v"),paste(uu[2],uu[1],sep = "v"))
    results<-decideTests(fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5)
    colnames(results)<-contrast.names[2]
}

if (opt$model == "Common Reference"){

    uu<-unique(targets2[,2])
    uuu1<-unique(targets3[,2])
    uuu2<-unique(targets3[,3]) #confirm reference column
    design <- modelMatrix(targets3, ref = uuu2)
    fit <- lmFit(MA.summarized, design)
    
    f <- factor(targets3[,2], levels=unique(uuu1)) #confirm non reference column
    
    # Create Contrast Model
    combos<-combn(levels(f),2) # generate matrix of pairwise group combinations for comparison
    ff <- factor(targets2[,2], levels=unique(uu))
    combos.names<-combn(levels(ff),2)
    contrasts<-c(paste(combos[1,],combos[2,],sep = "-"),paste(combos[2,],combos[1,],sep = "-")) # format combinations for limma:makeContrasts
    contrast.names <-c(paste(combos.names[1,],combos.names[2,],sep = "v"),paste(combos.names[2,],combos.names[1,],sep = "v")) # format combinations for output table file names
    
    
    cont.matrix <- makeContrasts(contrasts = contrasts,levels=design)
    contrast.fit <- contrasts.fit(fit, cont.matrix)
    contrast.fit <- eBayes(contrast.fit)
    results<-decideTests(contrast.fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5) # FDR .05
}
  
if (opt$model == "Separate Channels"){
    targets$t2 <- targetsA2C(targets$t2, channel.codes = c(1,2), channel.columns = list(Target=c("Cy3","Cy5")),grep = FALSE)
    targets$t3 <- targetsA2C(targets$t3, channel.codes = c(1,2), channel.columns = list(Target=c("Cy3","Cy5")),grep = FALSE)
    
    uuu <- unique(targets$t3$Target)
    uu <- unique(targets$t2$Target)
    fff <- factor(targets$t3$Target, levels=uuu)
    design <- model.matrix(~0+fff)
    colnames(design) <- uuu
    corfit <- intraspotCorrelation(MA.summarized, design)
    fit <- lmscFit(MA.summarized, design, correlation=corfit$consensus)
    if (opt$filter_nonexpressed){
      
      CutOff <- quantile(fit$Amean,probs=.33)
      
      hist_res <- graphics::hist(as.matrix(fit$coefficients), 100, col = "cornsilk", freq = FALSE, 
                                 main = "Probe Filtering Intensity Cutoff",
                                 border = "antiquewhite4",
                                 xlab = "Median intensities")
      
      abline(v = CutOff, col = "coral4", lwd = 2)
      
      keep <- fit$Amean > CutOff
      fit <- fit[keep,] # filter out probes below cutoff expression level
      annotation.subset <- annotation.subset[keep,]
      MA.summarized <- MA.summarized[keep,]
    }
    
    # Create Contrast Model
    combos<-combn(levels(fff),2) # generate matrix of pairwise group combinations for comparison
    uu <- unique(targets$t2$Target)
    ff <- factor(targets$t2$Target, levels=uu)
    combos.names<-combn(levels(ff),2)
    contrasts<-c(paste(combos[1,],combos[2,],sep = "-"),paste(combos[2,],combos[1,],sep = "-")) # format combinations for limma:makeContrasts
    contrast.names <-c(paste(combos.names[1,],combos.names[2,],sep = "v"),paste(combos.names[2,],combos.names[1,],sep = "v")) # format combinations for output table file names
    cont.matrix <- makeContrasts(contrasts = contrasts,levels=design)
    contrast.fit <- contrasts.fit(fit, cont.matrix)
    contrast.fit <- eBayes(contrast.fit, trend = TRUE)
    results<-decideTests(contrast.fit, method = "separate", adjust.method = "BH", p.value = 0.05, lfc = 0.5)
    
}

### Generate DGE Output Files

dir.create(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"), showWarnings = FALSE)
setwd(file.path(workdir,"Processed_Data",opt$glds,"02-Limma_DGE"))

########## Construct Output Data Tables
# output_table <- annotation.subset[,-c(1)]
# reduced_output_table <- annotation.subset[,-c(1)]
output_table <- annotation.subset
reduced_output_table <- annotation.subset
output_table <- cbind(output_table,as.data.frame(MA.summarized$A))
reduced_output_table <- cbind(reduced_output_table,as.data.frame(MA.summarized$A))


output_table$All.mean <- contrast.fit$Amean
reduced_output_table$All.mean <- contrast.fit$Amean

output_table$All.stdev <- contrast.fit$s2.post
reduced_output_table$All.stdev <- contrast.fit$s2.post

output_table$F.p.value <- contrast.fit$F.p.value
reduced_output_table$F.p.value <- contrast.fit$F.p.value

if (opt$model == "Separate Channels"){
  
  ########## Add Group Mean Values
  group_means<-as.data.frame(fit$coefficients)
  colnames(group_means)<-paste0("Group.Mean_",uu)
  output_table<-cbind(output_table,group_means)
  reduced_output_table<-cbind(reduced_output_table,group_means)
  rm(group_means)
  
  # add group stdev columns
  group_stdev<-as.data.frame(fit$stdev.unscaled * fit$coefficients)
  colnames(group_stdev)<-paste0("Group.Stdev_",uu)
  output_table<-cbind(output_table,group_stdev)
  reduced_output_table<-cbind(reduced_output_table,group_stdev)
  rm(group_stdev)  
  
}

if (opt$model == "Replicate Array"){
  uu<-uu[1]
  ########## Add Group Mean Values
  group_means<-as.data.frame(fit$coefficients)
  colnames(group_means)<-paste0("Group.Mean_",uu)
  output_table<-cbind(output_table,group_means)
  reduced_output_table<-cbind(reduced_output_table,group_means)
  rm(group_means)
  
  # add group stdev columns
  group_stdev<-as.data.frame(fit$stdev.unscaled * fit$coefficients)
  colnames(group_stdev)<-paste0("Group.Stdev_",uu)
  output_table<-cbind(output_table,group_stdev)
  reduced_output_table<-cbind(reduced_output_table,group_stdev)
  rm(group_stdev)  
  
}


if (opt$model == "Replicate Array"){
  # Contrast 1
  top <- topTable(fit, coef = 1, number = Inf, genelist = fit$genes$ID, adjust.method = "BH", sort.by = "none")
  table <- top[,c(2,5,6)] # Pull columns for Log2fc, P.value, Adj.p.value
  colnames(table)<- c("Log2fc","P.value","Adj.p.value")
  table.reduced <- table
  table$Updown <- sign(top$logFC)
  table$Sig.1 <- top$adj.P.Val<=0.1
  table$Sig.05 <- top$adj.P.Val<=0.05
  table$Log2_P.value <- log2(top$P.Value) # For volcano plot
  table$Log2_Adj.p.value <- log2(top$adj.P.Val) # For volcano plot
  colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[1],sep = "_")
  colnames(table)<-paste(colnames(table),contrast.names[1],sep = "_")
  output_table<-cbind(output_table,table)
  reduced_output_table<-cbind(reduced_output_table,table.reduced)
  
  # Repeat for contrast 2
  table <- top[,c(2,5,6)] # Pull columns for Log2fc, P.value, Adj.p.value
  colnames(table)<- c("Log2fc","P.value","Adj.p.value")
  table$Log2fc <- -1*table$Log2fc
  table.reduced <- table
  table$Updown <- sign(table$Log2fc)
  table$Sig.1 <- top$adj.P.Val<=0.1
  table$Sig.05 <- top$adj.P.Val<=0.05
  table$Log2_P.value <- log2(top$P.Value) # For volcano plot
  table$Log2_Adj.p.value <- log2(top$adj.P.Val) # For volcano plot
  
  colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[2],sep = "_")
  colnames(table)<-paste(colnames(table),contrast.names[2],sep = "_")
  output_table<-cbind(output_table,table)
  reduced_output_table<-cbind(reduced_output_table,table.reduced)
}

if (opt$model == "Separate Channels" || opt$model == "Common Reference"){
  
  # iterate through contrasts
  for (i in 1:length(contrasts)){
    top <- topTable(contrast.fit, coef = i, number = Inf, genelist = contrast.fit$genes$ID, adjust.method = "BH", sort.by = "none")
    table <- top[,c(2,5,6)] # Pull columns for Log2fc, P.value, Adj.p.value
    colnames(table)<- c("Log2fc","P.value","Adj.p.value")
    table.reduced <- table
    table$Updown <- sign(top$logFC)
    table$Sig.1 <- top$adj.P.Val<=0.1
    table$Sig.05 <- top$adj.P.Val<=0.05
    table$Log2_P.value <- log2(top$P.Value) # For volcano plot
    table$Log2_Adj.p.value <- log2(top$adj.P.Val) # For volcano plot
    
    colnames(table.reduced)<-paste(colnames(table.reduced),contrast.names[i],sep = "_")
    colnames(table)<-paste(colnames(table),contrast.names[i],sep = "_")
    output_table<-cbind(output_table,table)
    reduced_output_table<-cbind(reduced_output_table,table.reduced)
  }
  rm(i,top,table,table.reduced)
}

### Export DGE Output Data Tables


write.csv(reduced_output_table,"differential_expression.csv", row.names = FALSE)
write.csv(output_table,"visualization_output_table.csv", row.names = FALSE)


if (opt$model == "Replicate Array"){
  contrast.names<-c(paste(uu[1],uu[2],sep = "v"),paste(uu[2],uu[1],sep = "v"))
  write.csv(contrast.names,"contrasts.csv")
} else {
  contrast.output <- contrast.fit$contrasts
  row.names(contrast.output)<-uu
  colnames(contrast.output)<-contrast.names
  write.csv(contrast.output,"contrasts.csv")
}

### Export Metadata files
dir.create(file.path(workdir,"Processed_Data",opt$glds,"Metadata"), showWarnings = FALSE)
path<-file.path(workdir,"Processed_Data",opt$glds,"Metadata")
setwd(path)
file.copy(from = opt$isafile$datapath, to = file.path(path,opt$isafile$name),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE)
try(file.copy(from = opt$probefile$datapath, to = file.path(path,opt$probefile$name),overwrite = FALSE, recursive = FALSE, copy.mode = FALSE))

cat("All data files have been written to:  ",file.path(workdir,"Processed_Data"))