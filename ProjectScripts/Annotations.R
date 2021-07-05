BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")

AnnoLoc = "data/Annotations/"

AssemblyFilename = "gencode.v35.annotation.gff3.gz"
#AssemblyFilename = "gencode.v35lift37.annotation.gff3.gz"

txdbFilename  = paste0(paste0(strsplit(AssemblyFilename, "\\.")[[1]][c(2,4)], collapse = "."), "_DB")
AnnoFilename = paste0("annoFileCollapsed_", strsplit(txdbFilename, "\\.")[[1]][1], ".Rds")

#Get genome annotations

#annoFileCollapsed <- GetGenomeAnno_UCSC(genome = "hg19")
if(!file.exists(paste0(AnnoLoc, txdbFilename))){
  txdb <- makeTxDbFromGFF(paste0(AnnoLoc, AssemblyFilename))
  saveDb(txdb, paste0(AnnoLoc, txdbFilename))
} else {
  txdb <- loadDb(paste0(AnnoLoc, txdbFilename))
}

if(!file.exists(paste0(AnnoLoc, AnnoFilename))){
  annoFileCollapsed <- GetGenomeAnno_GENCODE(GTF_file = paste0(AnnoLoc, AssemblyFilename), txdb = txdb)
  saveRDS(annoFileCollapsed, paste0(AnnoLoc, AnnoFilename))
} else {
  annoFileCollapsed <- readRDS(paste0(AnnoLoc, AnnoFilename))
}
