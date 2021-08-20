#BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")
packageF("rvest")
AnnoLoc = "../Annotations/"

#AssemblyFilename = "gencode.v35.annotation.gff3.gz"
AssemblyFilename = "gencode.v35lift37.annotation.gff3.gz"
if(!file.exists(paste0(AnnoLoc, AssemblyFilename))){
  shortName = strsplit(AssemblyFilename, "\\.")[[1]][2]
  HTML <- read_html(paste0("https://www.gencodegenes.org/human/release_", gsub("^v", "", shortName), ".html"))
  URL <-  html_nodes(temp, "a") %>% .[grepl("GFF3", html_nodes(temp, "a") %>% html_text())] %>%
    as.character() %>% .[1]
  
  URL <- strsplit(URL, '"')[[1]][2]
  download.file(URL, destfile = paste0(AnnoLoc, AssemblyFilename))
}

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
