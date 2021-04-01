BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

source("ProjectScripts/ProjectFunctions.R")

GTF_file = paste0(AnnoLoc, AssemblyFilename)
geneNameFile = paste0(AnnoLoc, strsplit(AssemblyFilename, "\\.")[[1]][2], "_geneNames.Rds")

ResultsPath = paste0("Results_", Cohort, "_", strsplit(AssemblyFilename, "\\.")[[1]][2], "/RNAseq")

if(!ResultsPath %in% list.dirs(full.names = F, recursive = T)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

#Get transcript  annotations
if(!file.exists(geneNameFile)){
  geneNames <- rtracklayer::import.gff3(GTF_file) %>%
    data.frame() %>% filter(type == "transcript")  %>%
    select(seqnames, start, end, ID,
           gene_id, gene_name, gene_type, transcript_type)
  saveRDS(geneNames, geneNameFile)
} else {
  geneNames <- readRDS(geneNameFile)
}

#Get metadata
Metadata <- readRDS("meta/RNAseqMeta.Rds")

#Filter relevant samples
Metadata %<>% filter(Age > 15, Cohort != "Pool") %>% droplevels()
Metadata$Agef <- cut(Metadata$Age, breaks = 5, ordered_result = T)
Metadata %<>% mutate(Batch = paste0("Batch_", Batch)) 

Metadata$OrgRegion = factor("Cortex")
Metadata %<>% mutate(NeuExpRegion = OrgRegion,
                     Filename = RNAseq_id_ParkOme2,
                     Series_sample_id = RNAseq_id_ParkOme2,
                     Study = "Aging")


# Get counts
CountTxt <- readRDS("data/tximport_tx_hg19.Rds") %>%
  .$counts %>% data.frame() %>% 
  mutate(txt_ID_org = rownames(.))

CountTxt$txt_ID <- sapply(CountTxt$txt_ID_org, function(x){
  strsplit(x, "_")[[1]][1]
})

CountTxt$gene_id <- geneNames$gene_id[match(CountTxt$txt_ID, geneNames$ID)]

countMatrix <- CountTxt %>% 
  group_by(gene_id) %>% summarise(across(matches("SL"), sum)) %>% data.frame()

rownames(countMatrix) <- countMatrix$gene_id

countMatrix <- countMatrix[-1] %>% as.matrix()

# Match count matrix names to metadata names
countMatrix <- countMatrix[,Metadata$RNAseq_id_ParkOme2] 

#Remove genes with maximal count < 2 and mitochondrial genes
Max <- apply(countMatrix, 1, max)
mitoGenes <- geneNames %>% filter(seqnames == "chrM")

countMatrix <- countMatrix[Max > 10,]

countMito <- countMatrix[rownames(countMatrix) %in% mitoGenes$gene_id,]
countMitoSum <- countMito %>% apply(2, sum)
TotLibSize <- countMatrix  %>% apply(2, sum) 

MitoCountFiltered <- countMatrix[!rownames(countMatrix) %in% mitoGenes$gene_id,]
MitoFiltCountSum = apply(MitoCountFiltered, 2, sum)

Metadata$MitoCount <- countMitoSum[match(Metadata$RNAseq_id_ParkOme2, names(countMitoSum))]
Metadata$TotLibSize <- TotLibSize[match(Metadata$RNAseq_id_ParkOme2, names(TotLibSize))]

#Look at sample correlation after mitochndrial gene removal 
SampleCor <- cor(MitoCountFiltered)

annoCol = data.frame(Age = Metadata$Age,
                     Batch = Metadata$Batch,
                     Sex = Metadata$Sex,
                     Cohort = Metadata$Cohort,
                     MitoCount = Metadata$MitoCount,
                     LibSize = Metadata$TotLibSize,
                     row.names = Metadata$RNAseq_id_ParkOme2)

annoColors <- list(Cohort = c("Netherlands Brain Bank" = "dodgerblue4" , "Neuromics Tissue Bank" = "chocolate1"),
                   Sex = c(F = "indianred4", M = "cornflowerblue"))

annoColors <- sapply(names(annoColors), function(x){
  annoColors[[x]][names(annoColors[[x]]) %in% annoCol[[x]]]
},simplify = F)

annoColors$Age <- c("darkseagreen1", "darkorchid4")

pheatmap(SampleCor, angle_col = 90, na_col = "white",border_color = NA,
         color = colorRampPalette(c("darkblue", "gold2"))(999),
         annotation_col = annoCol,
         annotation_colors = annoColors)


#Look at the most hiighly expressed genes
TopFiveProportion <- sapply(names(MitoFiltCountSum), function(sbj){
  SubMatrix = data.frame(genes = rownames(MitoCountFiltered), Counts = MitoCountFiltered[,sbj])
  TopFive = SubMatrix %>% arrange(desc(Counts)) %>% head(5)
  TopFive %<>%  mutate(Proportion = Counts/MitoFiltCountSum[sbj])
  Genes <- geneNames[match(TopFive$genes, geneNames$gene_id),]  %>% select(gene_name, gene_type)
  Genes$Filename = sbj
  temp <- cbind(Genes, TopFive)
  names(temp)[names(temp) == "genes"] <- "ensemblID"
  temp
}, simplify = FALSE) %>% rbindlist()


TopFiveSum <- TopFiveProportion %>% group_by(Filename) %>%
  summarise(TotProp = sum(Proportion)) %>%
  data.frame %>% arrange(TotProp)


TopFiveGeneFreq <- TopFiveProportion %>% group_by(ensemblID) %>%
  summarise(n = n()) %>%
  data.frame

TopFiveGeneFreq <- merge(TopFiveGeneFreq, geneNames,
                         by.x = "ensemblID", by.y = "gene_id",
                         all.x = T, all.y = F, sort = F) 

TopFiveGeneFreq$ensemblID2 <- sapply(TopFiveGeneFreq$ensemblID,  function(x){
 strsplit(x, "\\.")[[1]][1]
})

TopFiveGeneFreq %<>% mutate(ensemblID2 = paste0(ensemblID2,
                                                " (", gene_name, ", ", n, ")"))
TopFiveGeneFreq %<>% arrange(desc(n))

TopFiveProportion$ensemblID2 <- TopFiveGeneFreq$ensemblID2[match(TopFiveProportion$ensemblID,
                                                                 TopFiveGeneFreq$ensemblID)]

#Order subjects based on library size after filtering of mtDNA genes
TopFiveProportion <- merge(TopFiveProportion, Metadata %>%
                             select(Filename, Batch, Cohort,
                                    DV200, RIN, MitoCount),
                           by = "Filename", sort = F)

TopFiveProportion$Filename <- factor(TopFiveProportion$Filename,
                                     levels = TopFiveSum$Filename)
TopFiveProportion$ensemblID2 <- factor(TopFiveProportion$ensemblID2,
                                       levels = unique(as.character(TopFiveGeneFreq$ensemblID2)))



TopFivePlot <- ggplot(TopFiveProportion %>% filter(Proportion > 0.006), aes(Filename, Proportion, fill = ensemblID2)) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  labs(y = "Proportion of reads", x = "Sample", title = "Top five genes with the highest read count") + 
  scale_fill_manual(values = c(MoviePalettes$MoonRiseKingdomColors[3:9], gray.colors(6)),
                    name = "GeneID (n)") +
  geom_bar(stat = "identity")

ggarrange(TopFivePlot,
          ggplot(TopFiveProportion %>% group_by(Filename) %>%
                                summarise(TopFiveProp = sum(Proportion),
                                          DV200 = mean(DV200),
                                          mtDNAcount = mean(MitoCount)) %>% data.frame(),
                              aes(DV200, TopFiveProp, color = mtDNAcount)) +
            geom_point(),
          nrow = 2, heights = c(1.5, 1))
ggsave(paste0(ResultsPath, "TopFiveGenesProp.pdf"),
       device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F)


#Get the common genes with the highest count in majority of the samples and rmove them from the count matrix
CommonTopGenes <- TopFiveGeneFreq[TopFiveGeneFreq$n > 0.5*ncol(MitoCountFiltered),] %>% filter(!duplicated(ensemblID))

CommonTopGenesSum <- apply(MitoCountFiltered[rownames(MitoCountFiltered )%in% CommonTopGenes$ensemblID,], 2, sum)

countMatrixFiltered <- MitoCountFiltered[!rownames(MitoCountFiltered) %in% as.character(CommonTopGenes$ensemblID),]

# Remove gene with 0 counts in  > 80% of the samples
ZeroCount <- apply(countMatrixFiltered, 1, function(x){
  sum(x==0)
})

countMatrixFiltered <- countMatrixFiltered[ZeroCount < 0.8*ncol(countMatrixFiltered),]

#Create log2 CPM matrix after removal of mitochondria-encoded genes
cpmMatrixFiltered <- Count2CPM(countMatrixFiltered) %>% data.frame()
cpmMatrixFiltered <- apply(cpmMatrixFiltered, c(1,2), function(x) log2(x+1)) %>% data.frame()
cpmMatrixFiltered <- cbind(rownames(countMatrixFiltered), cpmMatrixFiltered)
colnames(cpmMatrixFiltered)[1] <- "genes"


#Add gene symbols
GeneSymbolAll <- data.frame(GeneSymbol = geneNames$gene_name[match(rownames(cpmMatrixFiltered), geneNames$gene_id)],
                            Probe = rownames(cpmMatrixFiltered),
                            ensemblID = rownames(cpmMatrixFiltered))


ExpDataCPM <- cbind(GeneSymbolAll, cpmMatrixFiltered[-1])

#Look at sample correlation a
SampleCor2 <- cor(cpmMatrixFiltered %>% select(matches("SL")))
diag(SampleCor2) <- NA

annoCol = data.frame(Age = Metadata$Age,
                     Batch = Metadata$Batch,
                     Sex = Metadata$Sex,
                     Cohort = Metadata$Cohort,
                     MitoCount = Metadata$MitoCount,
                     LibSize = Metadata$TotLibSize,
                     row.names = Metadata$RNAseq_id_ParkOme2)

annoColors <- list(Cohort = c("Netherlands Brain Bank" = "dodgerblue4" , "Neuromics Tissue Bank" = "chocolate1"),
                   Sex = c(F = "indianred4", M = "cornflowerblue"))

annoColors <- sapply(names(annoColors), function(x){
  annoColors[[x]][names(annoColors[[x]]) %in% annoCol[[x]]]
},simplify = F)

annoColors$Age <- c("darkseagreen1", "darkorchid4")


pheatmap(SampleCor2, angle_col = 90, na_col = "white",border_color = NA,
         color = colorRampPalette(c("darkblue", "gold2"))(999),
         annotation_col = annoCol,
         annotation_colors = annoColors, filename = paste0(ResultsPath, "SampleCor.pdf"),
         width = 12, height = 10)
closeDev()



studyFinal <- PreProccessRNAseq(Metadata = Metadata, expData = ExpDataCPM,
                                SexCol = "Sex", Combat = FALSE, resultsPath = ResultsPath)
studyFinal$Metadata$PMI[is.na(studyFinal$Metadata$PMI)] <- median(studyFinal$Metadata$PMI,
                                                                  na.rm = T)
studyFinal$Metadata$PMI_binned <- cut(studyFinal$Metadata$PMI,
                                      breaks=c(seq(from=0, to=48, by=6),
                                               max(studyFinal$Metadata$PMI)),
                                      include.lowest=TRUE, labels=1:9) %>% as.integer

studyFinal$countMatrix <- countMatrixFiltered[rownames(countMatrixFiltered) %in% studyFinal$ExpHigh$ensemblID,]
studyFinal$countMatrix <- apply(studyFinal$countMatrix, c(1,2), as.integer)

##### Add cell composition estimates ########
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/Cell_type_PCA.R?raw=T")

region = studyFinal$Metadata$NeuExpRegion %>% unique
CellType_genes <- GetMarkers(region)

#Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]

#Bootstrap with replacement the samples (90% of the samples)
SampleNames <- as.character(studyFinal$Metadata$Filename)

PCAresults <- sapply(paste0("Boot_", 1:100), function(boot){
  BootSamples <- sample(SampleNames, 0.9*length(SampleNames), replace = F)
  dataSub <- studyFinal$ExpHigh %>% select(c("GeneSymbol", BootSamples))
  PCA_genes_All_based(dataset_id="Aging",
                      dataset=dataSub,
                      CellType_genes=CellType_genes,
                      contName = "SL",SampleReg = "SL",
                      NoiseThershold = studyFinal$NoiseThreshold)
}, simplify = F)

PCA_resultsMean <- sapply(names(PCAresults[[1]]$modified), function(celltype){
  temp <- data.frame(CommonName = names(PCAresults[[1]]$modified[[celltype]]$x[,1]),
                     Rot = PCAresults[[1]]$modified[[celltype]]$x[,1])
  for(i in 2:length(PCAresults)){
    temp <- merge(temp, data.frame(CommonName = names(PCAresults[[i]]$modified[[celltype]]$x[,1]),
                                   Rot = PCAresults[[i]]$modified[[celltype]]$x[,1]), by = "CommonName", all = TRUE)
  }
  names(temp)[2:ncol(temp)] <- paste0("Rot", c(1:c(ncol(temp)-1)))
  temp$MeanRot <- rowMeans(temp[-1], na.rm = T)
  temp
}, simplify=FALSE)


#Add estimation to Metadata 
AllEstimates <- lapply(PCA_resultsMean, function(x){
  x$MeanRot
}) %>% do.call(cbind, .) %>% data.frame()
AllEstimates$CommonName <- PCA_resultsMean[[1]]$CommonName

studyFinal$Metadata <- merge(studyFinal$Metadata,
                             AllEstimates, by = "CommonName",
                             sort = F)


##### Add estimates of synapses based on GO annotations
packageF("fgsea")
packageF("tibble")
PathwaysList <- list(Hallmark = gmtPathways("data/Enrichment/h.all.v7.2.symbols.xls"),
                     GObp = gmtPathways("data/Enrichment/c5.go.bp.v7.2.symbols.xls"),
                     GOmf = gmtPathways("data/Enrichment/c5.go.mf.v7.2.symbols.xls"),
                     GOcc = gmtPathways("data/Enrichment/c5.go.cc.v7.2.symbols.xls"),
                     KEGG = gmtPathways("data/Enrichment/c2.cp.kegg.v7.2.symbols.xls"))


PCApresynapse <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                       GOterm = PathwaysList$GOcc$GO_PRESYNAPSE)

PCApostsynapse <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                        GOterm = PathwaysList$GOcc$GO_POSTSYNAPSE)

PrePostIntersect <- intersect(PCApresynapse$GeneIn,
                              PCApostsynapse$GeneIn)

PreSynUnique <- PCApresynapse$GeneIn[!PCApresynapse$GeneIn %in% PrePostIntersect]

PCApresynapseUnique <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                             GOterm = PreSynUnique)

PostSynUnique <- PCApostsynapse$GeneIn[!PCApostsynapse$GeneIn %in% PrePostIntersect]

PCApostynapseUnique <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                             GOterm = PostSynUnique)

studyFinal$Metadata$PostSynapse_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                        MGP = PCApostsynapse$PCA$x[,1])
studyFinal$Metadata$PreSynapse_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                       MGP = PCApresynapse$PCA$x[,1])
studyFinal$Metadata$PostSynapseUnique_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                              MGP = PCApostynapseUnique$PCA$x[,1])
studyFinal$Metadata$PreSynapseUnique_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                             MGP = PCApresynapseUnique$PCA$x[,1])


# Check the main sources of variation
PCAall <- prcomp(t(studyFinal$ExpAll %>% select(matches("SL"))), scale = T)
studyFinal$Metadata <- merge(studyFinal$Metadata, PCAall$x[,c(1:5)],
                             by.x = "RNAseq_id_ParkOme2", by.y = "row.names", sort = F)



lm(PC1~Age + Sex + DV200 + Oligo_Genes +
     Endothelial_Genes + Microglia_Genes +
     GabaPV_Genes + Astrocyte_Genes + Pyramidal_Genes +
     PostSynapseUnique_Genes + PreSynapseUnique_Genes +
     Cohort, data = studyFinal$Metadata) %>% summary()

lm(Age~ Sex + DV200 + Oligo_Genes +
     Endothelial_Genes + Microglia_Genes +
     GabaPV_Genes + Astrocyte_Genes + Pyramidal_Genes +
     PostSynapseUnique_Genes + PreSynapseUnique_Genes +
     Cohort, data = studyFinal$Metadata) %>% summary()

#Remove columns for cell types that were not estimated
studyFinal$Metadata <- studyFinal$Metadata[,!apply(studyFinal$Metadata, 2, function(x){
  sum(is.na(x)) == nrow(studyFinal$Metadata)
})]

AgeChanges <- sapply(names(studyFinal$Metadata)[grepl("_Genes",
                                               names(studyFinal$Metadata))],
                     function(celltype){
                       lm(as.formula(paste0(celltype, "~Age + Sex + Cohort + DV200")),
                          data = studyFinal$Metadata) %>% summary %>% .$coef
                     }, simplify = F)

AgeChangesDF <- sapply(names(AgeChanges), function(CellType){
  data <- AgeChanges[[CellType]] %>% data.frame()
  names(data)[4] <- "pValue"
  data[2,c(1,4)] %>% data.frame() %>% mutate(CellType = CellType)
}, simplify = F) %>% rbindlist()

AgeChanges2 <- sapply(names(studyFinal$Metadata)[grepl("_Genes",
                                                      names(studyFinal$Metadata))],
                     function(celltype){
                       lm(as.formula(paste0(celltype, "~Age + Sex + Cohort + DV200 + Microglia_Genes")),
                          data = studyFinal$Metadata) %>% summary %>% .$coef
                     }, simplify = F)

AgeChangesDF2 <- sapply(names(AgeChanges2), function(CellType){
  data <- AgeChanges2[[CellType]] %>% data.frame()
  names(data)[4] <- "pValue"
  data[2,c(1,4)] %>% data.frame() %>% mutate(CellType = CellType)
}, simplify = F) %>% rbindlist()


write.table(AgeChangesDF, paste0(ResultsPath, "CellTypeAgeRNA.tsv"), sep = "\t",
            row.names = F, col.names = T)

rm(countMatrix, ExpDataCPM, GeneSymbolAll, CountTxt, countMatrix,
   mitoGenes, countMito, countMatrixFiltered, CommonTopGenes,
   MitoCountFiltered)


# #Run DE analysis
Model = as.formula(~Agef + Sex + DV200 + Cohort + Oligo_Genes + Microglia_Genes)

DESeqOut <- DESeq2runRNA(data =  studyFinal$countMatrix, Meta = studyFinal$Metadata, model = Model)

DESeqResults <- GetDESeq2ResultsRNA(DESeqOut, coef = "Agef.L")

Ranks = deframe(DESeqResults %>% filter(!duplicated(gene_name), !is.na(padj)) %>% select(gene_name, stat))

fgseaResultsRNAseq <- lapply(PathwaysList, function(PathType){
  fgseaMultilevel(pathways=PathType, stats=Ranks, nPermSimple = 1000)
})


Model2 = as.formula(~Agef + Sex + DV200 + Cohort + Oligo_Genes + Microglia_Genes + Pyramidal_Genes)

DESeqOut2 <- DESeq2RUN(data =  studyFinal$countMatrix, Meta = studyFinal$Metadata, model = Model2)

DESeqResults2 <- GetDESeq2Results(DESeqOut2, coef = "Agef.L")

Ranks2 = deframe(DESeqResults2 %>% filter(!duplicated(gene_name), !is.na(padj)) %>% select(gene_name, stat))

fgseaResultsRNAseq2 <- lapply(PathwaysList, function(PathType){
  fgseaMultilevel(pathways=PathType, stats=Ranks2, nPermSimple = 1000)
})

ChIPResults <- readRDS(paste0(DiscoveryResults, "/DESegResultsAge.L_FullAll.Rds"))
ChIPResultsOut <- readRDS(paste0(DiscoveryResults, "/DESeqOutAll_Full.Rds"))

#Selecting the most significant peak to represent a gene
temp <- merge(DESeqResults %>% select(gene_name, stat, pvalue, padj, gene_type, EnsemblID),
              ChIPResults %>% arrange(padj) %>%
                filter(!duplicated(gene_id)) %>%
                select(gene_id, stat, pvalue, padj),
              by.x = "EnsemblID", by.y = "gene_id", suffixes = c("_RNA", "_ChIP")) %>%
  mutate(DeltaStat = stat_ChIP - stat_RNA,
         SameDirect = stat_ChIP + stat_RNA) %>%
  filter(!duplicated(gene_name))

ggplot(temp %>% filter(gene_type == "protein_coding"), aes(stat_RNA, stat_ChIP)) +
  theme_minimal() +
  geom_density_2d_filled() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red")

#Lookinng only at promoters
temp2 <- merge(DESeqResults %>% select(gene_name, stat, pvalue, padj, gene_type, EnsemblID),
              ChIPResults %>% arrange(padj) %>% filter(region_type == "Promoters") %>%
                filter(!duplicated(gene_id)) %>%
                select(gene_id, stat, pvalue, padj),
              by.x = "EnsemblID", by.y = "gene_id", suffixes = c("_RNA", "_ChIP")) %>%
  mutate(DeltaStat = stat_ChIP - stat_RNA,
         SameDirect = stat_ChIP + stat_RNA) %>%
  filter(!duplicated(gene_name))

ggplot(temp2 %>% filter(gene_type == "protein_coding"), aes(stat_RNA, stat_ChIP)) +
  theme_minimal() +
  geom_density_2d_filled() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red")




Ranks = deframe(temp %>% filter(gene_type == "protein_coding") %>% select(gene_name, DeltaStat))

fgseaResultsMostSignif <- lapply(PathwaysList, function(PathType){
  fgseaMultilevel(pathways=PathType, stats=Ranks, nPermSimple = 10000)
})

Ranks2 = deframe(temp %>% filter(gene_type == "protein_coding") %>% select(gene_name, SameDirect))

fgseaResultsMostSignifSame <- lapply(PathwaysList, function(PathType){
  fgseaMultilevel(pathways=PathType, stats=Ranks2, nPermSimple = 10000)
})


packageF("rGREAT")

GetChIPenrich <- function(ChIPdata) {
  HyperBed <- ChIPdata %>% filter(!duplicated(PeakName), padj < 0.05, log2FoldChange > 0) %>%
    select(Peak.CHR, Peak.START, Peak.END, PeakName)
  
  HypoBed <- ChIPdata %>% filter(!duplicated(PeakName), padj < 0.05, log2FoldChange < 0) %>%
    select(Peak.CHR, Peak.START, Peak.END, PeakName)
  
  BackgroundRegionsBed <- ChIPdata %>% filter(!duplicated(PeakName),
                                              !is.na(padj)) %>%
    select(Peak.CHR, Peak.START, Peak.END, PeakName)
  
  HyperRegions <- submitGreatJob(HyperBed, bg = BackgroundRegionsBed, species = "hg19")
  GreatOut <- getEnrichmentTables(HyperRegions)
  GreatOut2 <- plotRegionGeneAssociationGraphs(HyperRegions, type = 1)
  
  
  HypoRegions <- submitGreatJob(HypoBed, bg = BackgroundRegionsBed, species = "hg19")
  GreatOutHypo <- getEnrichmentTables(HypoRegions)
  GreatOutHypo2 <- plotRegionGeneAssociationGraphs(HypoRegions, type = 1)
  return(list(HyperEnrich = GreatOut, HypoEnrich = GreatOutHypo,
              HyperAnno = GreatOut2, HypoAnno = GreatOutHypo2))
}


DiscoveryAgingnrich <- GetChIPenrich(ChIPResults)


save.image(paste0(ResultsPath,"RNAseq.RData"))
save(studyFinal, file = paste0(ResultsPath, "studyFinal.Rda"))
save(PCAresults, file = paste0(ResultsPath, "PCAresults.Rda"))




