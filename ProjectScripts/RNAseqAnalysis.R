#BiocManager::install("devtools")
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

annoColors <- list(Cohort = c("Netherlands Brain Bank" = MoviePalettes$BugsLife[6] ,
                              "Neuromics Tissue Bank" = MoviePalettes$BugsLife[4]),
                   Sex = c(F = "indianred4", M = "cornflowerblue"),
                   Age = c("grey90", MoviePalettes$BugsLife[8]),
                   DV200 = c("black", "orange"))

annoColors <- sapply(names(annoColors), function(x){
  annoColors[[x]][names(annoColors[[x]]) %in% annoCol[[x]]]
},simplify = F)


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
                     DV200 = Metadata$DV200,
                     Batch = Metadata$Batch,
                     Sex = Metadata$Sex,
                     Cohort = Metadata$Cohort,
                     MitoCount = Metadata$MitoCount,
                     LibSize = Metadata$TotLibSize,
                     row.names = Metadata$RNAseq_id_ParkOme2)


pheatmap(SampleCor2, angle_col = 90, na_col = "white",border_color = NA, clustering_method = "ward.D2",
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


#Get Human microglia expression from Galatro et al. 2017
MicrogliaAging <- read.table("data/MicrogliaHumanAgingGalatro.txt", header = T, sep = "\t")
MicogliaAll <- read.table("data/MicrogliaVsWholeBrainGalatro.txt", header = T, sep = "\t")

MicrogliaSpecificHuman <- MicogliaAll %>% filter(gliaVSbrain_logFC > 3, adj.P.Val < 0.05) %>% .$GeneSymbol

OverlapMicrogliaAll <- intersect(MicrogliaSpecificHuman, rownames(PCAresults$Boot_1$All$Microglia_Genes$rotation))
OverlapAgeUP <- intersect((MicrogliaAging %>% filter(logFC > 0) %>% .$GeneSymbol), rownames(PCAresults$Boot_1$All$Microglia_Genes$rotation))
OverlapAgeDown <- intersect((MicrogliaAging %>% filter(logFC < 0) %>% .$GeneSymbol), rownames(PCAresults$Boot_1$All$Microglia_Genes$rotation))

#Add estimation to Metadata 
AllEstimates <- lapply(PCA_resultsMean, function(x){
  x$MeanRot
}) %>% do.call(cbind, .) %>% data.frame()
AllEstimates$CommonName <- PCA_resultsMean[[1]]$CommonName

studyFinal$Metadata <- merge(studyFinal$Metadata,
                             AllEstimates, by = "CommonName",
                             sort = F)


##### Add estimates of synapses based on GO annotations
PathwaysList <- list(Hallmark = gmtPathways("data/Enrichment/h.all.v7.2.symbols.xls"),
                     GObp = gmtPathways("data/Enrichment/c5.go.bp.v7.2.symbols.xls"),
                     GOmf = gmtPathways("data/Enrichment/c5.go.mf.v7.2.symbols.xls"),
                     GOcc = gmtPathways("data/Enrichment/c5.go.cc.v7.2.symbols.xls"),
                     KEGG = gmtPathways("data/Enrichment/c2.cp.kegg.v7.2.symbols.xls"))

PathwaysList2 <- sapply(names(PathwaysList), function(GeneSetType){
  GeneSet <- PathwaysList[[GeneSetType]]
  names(GeneSet) <- sapply(names(GeneSet), function(Ptwy){
    gsub("^GO_|^KEGG_|^HALL.*_", "", Ptwy) %>% gsub(" ", "_", .)
  })
  names(GeneSet) <- paste0(GeneSetType, "__", names(GeneSet))
  GeneSet
}, simplify = F)

PathwaysCombined <- c(PathwaysList2$GObp,
                      PathwaysList2$GOmf,
                      PathwaysList2$GOcc,
                      PathwaysList2$MitoCarta,
                      PathwaysList2$KEGG)

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
     Cohort,
   data = studyFinal$Metadata) %>% summary()

lm(PC1~Microglia_Genes + PreSynapseUnique_Genes +
     PostSynapseUnique_Genes,
   data = studyFinal$Metadata) %>% summary()


lm(PC1~Astrocyte_Genes + PreSynapseUnique_Genes +
     PostSynapseUnique_Genes,
   data = studyFinal$Metadata) %>% summary()

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
write.table(AgeChangesDF, paste0(ResultsPath, "CellTypeAgeRNA.tsv"), sep = "\t",
            row.names = F, col.names = T)


AgeChangesCI <-  sapply(names(studyFinal$Metadata)[grepl("_Genes",
                                                         names(studyFinal$Metadata))],
                        function(celltype){
                          lmOut <- lm(as.formula(paste0(celltype, "~Age + Sex + Cohort + DV200")),
                                      data = studyFinal$Metadata)
                          temp <- data.frame(Coef = lmOut %>% summary %>% .$coef %>%  .[2,1], Type = "Demographic adjustment")
                          temp <- cbind(temp, t(confint(lmOut)[2,]) %>% data.frame())
                          names(temp)[3:4] <- c("Low", "High")
                          temp$CellType = gsub("Genes", "MGP", celltype)
                          temp
                        }, simplify = F) %>% rbindlist() %>% data.frame()

AgeChangesCI_corrected <-  sapply(names(studyFinal$Metadata)[grepl("_Genes",
                                                                   names(studyFinal$Metadata))],
                                  function(celltype){
                                    lmOut <- lm(as.formula(paste0(celltype, "~Age + Sex + Cohort +  + DV200 + Microglia_Genes")),
                                                data = studyFinal$Metadata)
                                    temp <- data.frame(Coef = lmOut %>% summary %>% .$coef %>%  .[2,1],  Type = "Demographic and microglia adjustment")
                                    temp <- cbind(temp, t(confint(lmOut)[2,]) %>% data.frame())
                                    names(temp)[3:4] <- c("Low", "High")
                                    temp$CellType = gsub("Genes", "MGP", celltype)
                                    temp
                                  }, simplify = F) %>% rbindlist() %>% data.frame()

AgeChangesCIcombined <- rbind(AgeChangesCI, AgeChangesCI_corrected)

AgeChangesCIcombined$CellType <- factor(AgeChangesCI$CellType, levels = c("Astrocyte_MGP", "Endothelial_MGP", "Microglia_MGP", "Microglia_activation_MGP",               
                                                                          "Microglia_deactivation_MGP", "Oligo_MGP", "OligoPrecursors_MGP",
                                                                          "GabaPV_MGP", "GabaRelnCalb_MGP", "GabaVIPReln_MGP", "Layer_6b_Pyra_MGP", 
                                                                          "Pyramidal_MGP", "PostSynapse_MGP", "PreSynapse_MGP",
                                                                          "PostSynapseUnique_MGP", "PreSynapseUnique_MGP"))

ggplot(AgeChangesCIcombined[!grepl("Synapse", AgeChangesCIcombined$CellType),] %>% droplevels() , aes(CellType, Coef)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "", y = "Coeficient (95%CI)", size = 16) +
  geom_point() +
  geom_errorbar(aes(ymin = Low, ymax = High)) +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  facet_wrap(~Type, nrow = 2)
ggsave(paste0(ResultsPath, "MGPoutput.pdf"), device = "pdf", width = 6, height = 5, dpi = 300, useDingbats = F)

rm(countMatrix, ExpDataCPM, GeneSymbolAll, CountTxt, countMatrix,
   mitoGenes, countMito, countMatrixFiltered, CommonTopGenes,
   MitoCountFiltered)


# Run DE analysis
Model = as.formula(~Agef + Sex + DV200 + Cohort + Oligo_Genes + Microglia_Genes)

DESeqOut <- DESeq2runRNA(data =  studyFinal$countMatrix, Meta = studyFinal$Metadata, model = Model)

DESeqResults <- GetDESeq2ResultsRNA(DESeqOut, coef = "Agef.L")
EnrichmentOut <- RunEnrich(DESeqResults, method = "ErmineJ")

# Repeat correcting for Microglia and Pyramidal genes
Model2 = as.formula(~Agef + Sex + DV200 + Cohort + Oligo_Genes + Microglia_Genes + Layer_6b_Pyra_Genes)

DESeqOut2 <- DESeq2runRNA(data =  studyFinal$countMatrix, Meta = studyFinal$Metadata, model = Model2)

DESeqResults2 <- GetDESeq2ResultsRNA(DESeqOut2, coef = "Agef.L")
EnrichmentOut2 <- RunEnrich(DESeqResults2, method = "ErmineJ")

GOccTop <- rbind(EnrichmentOut$ErmineJ$EnrichUp$results %>% arrange(CorrectedPvalue) %>%
                   filter(CorrectedPvalue < 0.05, GeneSet == "GOcc") %>% head(10) %>%
                   mutate(Mod = "Oligo and Microglia"),
                 EnrichmentOut$ErmineJ$EnrichDown$results %>% arrange(CorrectedPvalue) %>%
                   filter(CorrectedPvalue < 0.05, GeneSet == "GOcc") %>% head(10) %>%
                   mutate(Mod = "Oligo and Microglia"),
                 EnrichmentOut2$ErmineJ$EnrichUp$results %>% arrange(CorrectedPvalue) %>%
                   filter(CorrectedPvalue < 0.05, GeneSet == "GOcc") %>% head(10) %>%
                   mutate(Mod = "Oligo, Microglia and Pyramidal"),
                 EnrichmentOut2$ErmineJ$EnrichDown$results %>% arrange(CorrectedPvalue) %>%
                   filter(CorrectedPvalue < 0.05, GeneSet == "GOcc") %>% head(10) %>%
                   mutate(Mod = "Oligo, Microglia and Pyramidal")) 

GOccTop2 <- rbind(EnrichmentOut$ErmineJ$EnrichUp$results %>%
                    filter(CorrectedPvalue < 0.05, GeneSet == "GOcc", ID %in% GOccTop$ID) %>%
                    mutate(Mod = "Oligo and Microglia"),
                  EnrichmentOut$ErmineJ$EnrichDown$results %>%
                    filter(CorrectedPvalue < 0.05, GeneSet == "GOcc", ID %in% GOccTop$ID) %>%
                    mutate(Mod = "Oligo and Microglia"),
                  EnrichmentOut2$ErmineJ$EnrichUp$results %>%
                    filter(CorrectedPvalue < 0.05, GeneSet == "GOcc", ID %in% GOccTop$ID) %>%
                    mutate(Mod = "Oligo, Microglia and Pyramidal"),
                  EnrichmentOut2$ErmineJ$EnrichDown$results %>%
                    filter(CorrectedPvalue < 0.05, GeneSet == "GOcc", ID %in% GOccTop$ID) %>%
                    mutate(Mod = "Oligo, Microglia and Pyramidal")) %>% select(-GeneMembers) %>% data.frame()

GOccTopDispaly <- GOccTop2
GOccTopDispaly$RawScore <- apply(GOccTopDispaly %>% select(RawScore, Direction), 1, function(x){
  if(x[2] == "Down"){
    -1*as.numeric(x[1])
  } else{
    as.numeric(x[1])
  }
})

GOccTopDispaly %<>% arrange(RawScore)

GOccTopDispaly$ID <- factor(GOccTopDispaly$ID, levels = unique(GOccTopDispaly$ID))

ggplot(GOccTopDispaly, aes(ID, RawScore, fill = Direction)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  labs(x = "") +
  scale_fill_manual(values = MoviePalettes$SpiritedAway[c(3, 7)]) +
  coord_flip() +
  facet_wrap(~Mod)
ggsave(paste0(ResultsPath, "RNAseqEnrichTopCC.pdf"), device = "pdf", width = 10, height = 8)


#Run the analysis without microglia or pyramidal correction correction
Model3 = as.formula(~Agef + Sex + DV200 + Cohort + Oligo_Genes)

DESeqOut3 <- DESeq2runRNA(data =  studyFinal$countMatrix, Meta = studyFinal$Metadata, model = Model3)

DESeqResults3 <- GetDESeq2ResultsRNA(DESeqOut3, coef = "Agef.L")
EnrichmentOut3 <- RunEnrich(DESeqResults3, method = "ErmineJ")


#Combine with the ChIPseq results
ChIPResults <- readRDS(paste0(DiscoveryResults, "/DESegResultsAge.L_FullAll.Rds"))
ChIPResultsOut <- readRDS(paste0(DiscoveryResults, "/DESeqOutAll_Full.Rds"))
ChIPmeta <- colData(ChIPResultsOut) %>% data.frame
studyFinal$Metadata %<>% mutate(ChIPsampleID = paste0("X", ChIPseq_id))


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

temp2 <- merge(DESeqResults2 %>% select(gene_name, stat, pvalue, padj, gene_type, EnsemblID),
              ChIPResults %>% arrange(padj) %>%
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


#Repeat with the model without microglia correction
temp3 <- merge(DESeqResults3 %>% select(gene_name, stat, pvalue, padj, gene_type, EnsemblID),
              ChIPResults %>% arrange(padj) %>% filter(region_type == "Promoters") %>%
                select(gene_id, stat, pvalue, padj),
              by.x = "EnsemblID", by.y = "gene_id", suffixes = c("_RNA", "_ChIP")) %>%
  mutate(DeltaStat = stat_ChIP - stat_RNA,
         SameDirect = stat_ChIP + stat_RNA) %>%
  filter(!duplicated(gene_name))

ggplot(temp3 %>% filter(gene_type == "protein_coding"), aes(stat_RNA, stat_ChIP)) +
  theme_minimal() +
  geom_density_2d_filled() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red")

temp3 %>% filter(gene_name %in% PathwaysList$KEGG$KEGG_OXIDATIVE_PHOSPHORYLATION) %>%
  ggplot(aes(stat_RNA, stat_ChIP)) +
  theme_minimal() +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red")


temp3 %>% filter(gene_name %in% rownames(PCAresults$Boot_1$All$Microglia_Genes$rotation)) %>%
  ggplot(aes(stat_RNA, stat_ChIP)) +
  theme_minimal() +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red")

temp3 %>% filter(gene_name %in% rownames(PCAresults$Boot_1$All$Pyramidal_Genes$rotation)) %>%
  ggplot(aes(stat_RNA, stat_ChIP)) +
  theme_minimal() +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red")

#Look at the enrichment 
EnrichmentDelta <-  RunEnrich(temp %>% filter(gene_type == "protein_coding"),Filter = F, 
                              method = "ErmineJ", ScoreCol = "DeltaStat")
EnrichmentSum <-  RunEnrich(temp %>% filter(gene_type == "protein_coding"),Filter = F, 
                              method = "ErmineJ", ScoreCol = "SameDirect")

EnrichmentDelta2 <-  RunEnrich(temp2 %>% filter(gene_type == "protein_coding"),Filter = F, 
                              method = "ErmineJ", ScoreCol = "DeltaStat")
EnrichmentSu2m <-  RunEnrich(temp2 %>% filter(gene_type == "protein_coding"),Filter = F, 
                            method = "ErmineJ", ScoreCol = "SameDirect")


EnrichmentDelta3 <-  RunEnrich(temp3 %>% filter(gene_type == "protein_coding"),Filter = F, 
                              method = "ErmineJ", ScoreCol = "DeltaStat")
EnrichmentSum3 <-  RunEnrich(temp3 %>% filter(gene_type == "protein_coding"),Filter = F, 
                            method = "ErmineJ", ScoreCol = "SameDirect")


# CreateAdjCovar <- function(dds){
#   temp <- attr(dds, "modelMatrix") %>%
#     data.frame() %>% select(-matches("Age|Interc"))
#   DF <- data.frame(Cov = names(temp))
#   DF$adjType <- apply(temp, 2, function(Cov){
#     if(sum(as.integer(Cov) != Cov) == 0){
#       "base"
#     } else {
#       "mean"
#     }
#   })
#   return(DF)
# }
# 
# AdjCountsWrap <- function(dds, genes = NULL){
#   if(is.null(genes)){
#     genes = rownames(dds)
#   }
#   
#   adjCov <- CreateAdjCovar(dds)
# 
#   AdjRNA <- sapply(genes, function(gene){
#     temp <- GetAdjCountDESeq(dds = dds, gene, adjCov = adjCov)
#   })
#   rownames(AdjRNA) <- colnames(dds)
#   return(list(CovarDF = adjCov, AdjValDF = t(AdjRNA) %>% data.frame()))
# }
# 
# AdjCovRNA <- AdjCountsWrap(DESeqOut)
# saveRDS(AdjCovRNA, paste0(ResultsPath, "AgingAdjCovRNA.Rds"))
# 
# SubChIP <- ChIPResults %>%
#   filter(region_type == "Promoters",
#          !is.na(symbol),
#          gene_id %in% rownames(AdjCovRNA$AdjValDF), !is.na(padj)) %>%
#   filter(!duplicated(Peak_Gene))
# 
# AdjCovChIP <- AdjCountsWrap(ChIPResultsOut, genes = SubChIP %>%
#                               filter(!duplicated(PeakName)) %>% .$PeakName)
# saveRDS(AdjCovChIP, paste0(ResultsPath, "AgingAdjCovChIPseq.Rds"))
# 
# SlidingAge <- list()
# Start = 1
# End = 10
# 
# while(End <= 40){
#   ageStart = sort(studyFinal$Metadata$Age)[Start]
#   ageEnd = sort(studyFinal$Metadata$Age)[End]
#   SlidingAge[[paste0("Age",ageStart, "_", ageEnd)]] <- studyFinal$Metadata %>%
#     arrange(Age) %>% select(Age, Biobank_ID,
#                             RNAseq_id_ParkOme2, ChIPseq_id) %>%
#     .[Start:End,] %>% mutate(ChIPseq_id = paste0("X", ChIPseq_id))
#   
#   Start = Start + 5
#   End = End + 5
# }
# 
# SampleGroups = SlidingAge
# GenePairs <- rownames(AdjCovRNA$AdjValDF)[rownames(AdjCovRNA$AdjValDF) %in% SubChIP$gene_id]
# 
# # parallelize depending on the operating system
# packageF("doParallel")
# tempFun <- function(x){
#   library(magrittr)
#   library(dplyr)
#   library(data.table)
#   
#   AgeSamples = SampleGroups[[x]]
#   sapply(GenePairs, function(GeneID){
#     
#     subData = SubChIP[SubChIP$gene_id == GeneID,]
#     Peaks = subData %>% .$PeakName %>% unique
#     GeneType = subData %>% .$gene_type %>% unique
#     GeneSymbol = subData %>% .$symbol %>% unique
#     
#     PeakCor <- sapply(Peaks, function(PeakName){
#       RNAdata = AdjCovRNA$AdjValDF[rownames(AdjCovRNA$AdjValDF) == GeneID,] %>%
#         select(AgeSamples$RNAseq_id_ParkOme2) %>% unlist
#       ChIPdata = AdjCovChIP$AdjValDF[rownames(AdjCovChIP$AdjValDF) == PeakName,] %>%
#         select(AgeSamples$ChIPseq_id) %>% unlist
#       data.frame(GeneID = GeneID, PeakName = PeakName,
#                  GeneType = GeneType, GeneSymbol = GeneSymbol,
#                  Cor = cor(RNAdata, ChIPdata))
#     }, simplify = F) %>% rbindlist() %>% data.frame()
#     PeakCor %>% mutate(PairNum = nrow(.)) %>% arrange(desc(Cor))
#   }, simplify = F) %>% rbindlist %>% data.frame()
# }
# 
# if(.Platform$OS.type == "windows") {
#   
#   cl <- makeCluster(detectCores(), type = "PSOCK")  
#   registerDoParallel(cl)  
#   
#   clusterExport(cl, varlist = c("SampleGroups", "SubChIP",
#                                 "AdjCovRNA", "AdjCovChIP",
#                                 "GenePairs"))
#   
#   ChiP_RNAcorMax <- clusterApply(cl, x = names(SampleGroups), fun = tempFun)
#   names(ChiP_RNAcorMax) <- names(SampleGroups)
#   stopCluster(cl)
# } else {
#   ChiP_RNAcorMax <- mclapply(SampleGroups, function(AgeSamples){
#     sapply(GenePairs, function(GeneID){
#       subData = SubChIP[SubChIP$gene_id == GeneID,]
#       Peaks = subData %>% .$PeakName %>% unique
#       GeneType = subData %>% .$gene_type %>% unique
#       GeneSymbol = subData %>% .$gene_name %>% unique
#       PeakCor <- sapply(Peaks, function(PeakName){
#         RNAdata = AdjCovRNA$AdjValDF[rownames(AdjCovRNA$AdjValDF) == GeneID,] %>%
#           select(AgeSamples$RNAseq_id_ParkOme2) %>% unlist
#         ChIPdata = AdjCovChIP$AdjValDF[rownames(AdjCovChIP$AdjValDF) == PeakName,] %>%
#           select(AgeSamples$ChIPseq_id) %>% unlist
#         data.frame(GeneID = GeneID, PeakName = PeakName,
#                    GeneType = GeneType, GeneSymbol = GeneSymbol,
#                    Cor = cor(RNAdata, ChIPdata))
#       }, simplify = F) %>% rbindlist() %>% data.frame()
#       PeakCor %>% mutate(PairNum = nrow(.)) %>% arrange(desc(Cor))
#     }, simplify = F) %>% rbindlist %>% data.frame() 
#   }, mc.cores = detectCores())
# }
# 
# 
# ChiP_RNAcorMax <- sapply(names(ChiP_RNAcorMax), function(x){
#   ChiP_RNAcorMax[[x]] %>% mutate(AgeGroup = x,
#                                  GeneAge = paste0(GeneID, "_", AgeGroup),
#                                  GenePeak = paste0(GeneID, "_", PeakName))
# }, simplify = F)
# 
# saveRDS(ChiP_RNAcorMax, paste0(ResultsPath, "ChiP_RNAcorMax.Rds"))
# 
# CorAgeChange <- sapply(ChiP_RNAcorMax$Age17.9_55$GenePeak, function(GenePair){
#   Data = lapply(ChiP_RNAcorMax, function(AgeGroup){
#     AgeGroup %>% filter(GenePeak == GenePair)
#   }) %>% rbindlist() %>% data.frame()
#   Data$AgeGroupNum = 1:nrow(Data)
#   if(sum(!is.na(Data$Cor)) >= 4){
#     temp <- lm(Cor~AgeGroupNum, data = Data) %>%
#       summary() %>% .$coef %>% .[2,-2] %>% data.frame() %>% t()
#     colnames(temp) <- c("Est", "tStat", "pVal")
# 
#   } else {
#     temp = data.frame("Est" = NA, "tStat" = NA, "pVal" = NA)
#   }
# 
#   cbind(Data[1,] %>% select(GenePeak, GeneSymbol, PairNum), temp)
# 
# }, simplify = F) %>% rbindlist() %>% data.frame()
# 
# CorAgeChange$YoungCor = ChiP_RNAcorMax$Age17.9_55$Cor
# saveRDS(CorAgeChange, paste0(ResultsPath, "CorAgeChange.Rds"))
# 
# ggplot(CorAgeChange, aes(YoungCor, Est)) + geom_point() + geom_smooth()
# 
# 
# YoungHigh <- ChiP_RNAcorMax$Age17.9_55 %>% filter(Cor > 0.5)
# YoungLow <- ChiP_RNAcorMax$Age17.9_55 %>% filter(Cor < -0.5)
# 
# OldHigh <- ChiP_RNAcorMax$Age87_102 %>% filter(Cor > 0.4)
# OldLow <- ChiP_RNAcorMax$Age87_102 %>% filter(Cor < -0.4)
# 
# 
# 
# ggplot(rbindlist(ChiP_RNAcorMax) %>% data.frame() %>%
#   filter(GenePeak %in% YoungHigh$GenePeak)) + geom_boxplot(aes(AgeGroup, Cor))
# 
# ggplot(rbindlist(ChiP_RNAcorMax) %>% data.frame() %>%
#          filter(GenePeak %in% YoungLow$GenePeak)) + geom_boxplot(aes(AgeGroup, Cor))
# 
# ggplot(rbindlist(ChiP_RNAcorMax) %>% data.frame() %>%
#   filter(GenePeak %in% OldLow$GenePeak)) + geom_boxplot(aes(AgeGroup, Cor))
# 
# ggplot(rbindlist(ChiP_RNAcorMax) %>% data.frame() %>%
#          filter(GenePeak %in% OldHigh$GenePeak)) + geom_boxplot(aes(AgeGroup, Cor))
# 
# 
# #Create random correlations
# RandCorList <- list()
# for(i in c(1:10)){
#   SlidingAgeRand <- list()
#   Samples <- sample(studyFinal$Metadata$Biobank_ID, nrow(studyFinal$Metadata), replace = F)
#   Start = 1
#   End = 10
#   
#   while(End <= 40){
#     ageStart = studyFinal$Metadata[match(Samples, studyFinal$Metadata$Biobank_ID),]$Age[Start]
#     ageEnd = studyFinal$Metadata[match(Samples, studyFinal$Metadata$Biobank_ID),]$Age[End]
#     Data <- studyFinal$Metadata[match(Samples,
#                                       studyFinal$Metadata$Biobank_ID),] %>%
#       select(Age, Biobank_ID,
#              RNAseq_id_ParkOme2, ChIPseq_id) %>%
#       .[Start:End,] %>% mutate(ChIPseq_id = paste0("X", ChIPseq_id))
#     
#     SlidingAgeRand[[paste0("Age", min(Data$Age), "_", max(Data$Age))]] <- Data 
#     
#     Start = Start + 5
#     End = End + 5
#   }
#   
#   SampleGroups = SlidingAgeRand
#   
#   # parallelize depending on the operating system
#   if(.Platform$OS.type == "windows") {
#     
#     cl <- makeCluster(detectCores(), type = "PSOCK")  
#     registerDoParallel(cl)  
#     
#     clusterExport(cl, varlist = c("SampleGroups", "SubChIP",
#                                   "AdjCovRNA", "AdjCovChIP",
#                                   "GenePairs"))
#     
#     ChiP_RNAcorRand <- clusterApply(cl, x = names(SampleGroups), fun = tempFun)
#     names(ChiP_RNAcorRand) <- names(SampleGroups)
#     stopCluster(cl)
#   } else {
#     ChiP_RNAcorRand <- mclapply(SampleGroups, function(AgeSamples){
#       sapply(GenePairs, function(GeneID){
#         subData = SubChIP[SubChIP$gene_id == GeneID,]
#         Peaks = subData %>% .$PeakName %>% unique
#         GeneType = subData %>% .$gene_type %>% unique
#         GeneSymbol = subData %>% .$gene_name %>% unique
#         PeakCor <- sapply(Peaks, function(PeakName){
#           RNAdata = AdjCovRNA$AdjValDF[rownames(AdjCovRNA$AdjValDF) == GeneID,] %>%
#             select(AgeSamples$RNAseq_id_ParkOme2) %>% unlist
#           ChIPdata = AdjCovChIP$AdjValDF[rownames(AdjCovChIP$AdjValDF) == PeakName,] %>%
#             select(AgeSamples$ChIPseq_id) %>% unlist
#           data.frame(GeneID = GeneID, PeakName = PeakName,
#                      GeneType = GeneType, GeneSymbol = GeneSymbol,
#                      Cor = cor(RNAdata, ChIPdata))
#         }, simplify = F) %>% rbindlist() %>% data.frame()
#         PeakCor %>% mutate(PairNum = nrow(.)) %>% arrange(desc(Cor))
#       }, simplify = F) %>% rbindlist %>% data.frame() 
#     }, mc.cores = detectCores())
#   }
#   
#   
#   ChiP_RNAcorRand <- sapply(names(ChiP_RNAcorRand), function(x){
#     ChiP_RNAcorRand[[x]] %>% mutate(AgeGroup = x,
#                                     GeneAge = paste0(GeneID, "_", AgeGroup),
#                                     GenePeak = paste0(GeneID, "_", PeakName))
#   }, simplify = F)
#   
#   CorAgeChangeRand <- sapply(ChiP_RNAcorRand[[1]]$GenePeak, function(GenePair){
#     Data = lapply(ChiP_RNAcorRand, function(AgeGroup){
#       AgeGroup %>% filter(GenePeak == GenePair)
#     }) %>% rbindlist() %>% data.frame()
#     Data$AgeGroupNum = 1:nrow(Data)
#     if(sum(!is.na(Data$Cor)) >= 4){
#       temp <- lm(Cor~AgeGroupNum, data = Data) %>%
#         summary() %>% .$coef %>% .[2,-2] %>% data.frame() %>% t()
#       colnames(temp) <- c("Est", "tStat", "pVal")
#       
#     } else {
#       temp = data.frame("Est" = NA, "tStat" = NA, "pVal" = NA)
#     }
#     
#     cbind(Data[1,] %>% select(GenePeak, GeneSymbol, PairNum), temp)
#     
#   }, simplify = F) %>% rbindlist() %>% data.frame()
#   CorAgeChangeRand$YoungCor = ChiP_RNAcorRand[[1]]$Cor
#   RandCorList[[i]] <- list(SampleGroups = SlidingAgeRand,
#                            ChiP_RNAcorRand = ChiP_RNAcorRand,
#                            CorAgeChangeRand = CorAgeChangeRand)
# }
# 
# saveRDS(RandCorList, paste0(ResultsPath, "RandCorList.Rds"))


save.image(paste0(ResultsPath,"RNAseq.RData"))
save(studyFinal, file = paste0(ResultsPath, "studyFinal.Rda"))
save(PCAresults, file = paste0(ResultsPath, "PCAresults.Rda"))


lmAdjMicroglia <- lm(Microglia_Genes~Age+Sex+Cohort + DV200, data = studyFinal$Metadata)
studyFinal$Metadata$AdjMicroglia <- ModelAdj(lmAdjMicroglia,
                                             adj=data.frame(effect = c("Sex", "Cohort", "DV200"), adjValue=c(0, 0, 80)))

AdjMicrogliaPlot <- ggplot(studyFinal$Metadata, aes(Age, AdjMicroglia, color = Cohort)) +
  theme_minimal() +
  theme(legend.position = c(0.2,0.8)) +
  labs(y = "Microglia MGP (Adjusted)") +
  scale_color_manual(values = MoviePalettes$MadMaxDesert[c(1,7)]) +
  geom_smooth(color = "black", size = 0.5) +
  geom_point() 
  


lmModRiP <- lm(log(RiP_NormMeanRatioOrg)~Age + Sex + FinalBatch + NeuNall_MSP + Oligo_MSP, data = ChIPmeta)

ChIPmeta$AdjustedRiP <- ModelAdj(lmModRiP,
                                 adj=data.frame(effect = c("FinalBatch", "Sex",
                                                           "NeuNall_MSP", "Oligo_MSP"), adjValue=c(0, 0, 0.5, 0.5)))

AdjRiPPlot <- ggplot(ChIPmeta, aes(Age, AdjustedRiP, color = Cohort))+
  theme_minimal() +
  theme(legend.position = c(0.2,0.8)) +
  labs(y = "log(Normalized RiP) (Adjusted)") +
  scale_color_manual(values = MoviePalettes$MadMaxDesert[c(1,7)]) +
  geom_smooth(color = "black", size = 0.5) +
  geom_point(show.legend = F) 



studyFinal$Metadata$AdjustedRiP <- ChIPmeta$AdjustedRiP[match(studyFinal$Metadata$ChIPsampleID, ChIPmeta$SampleID)]

MicroRiPCorStat <- cor.test(~AdjMicroglia+AdjustedRiP, data = studyFinal$Metadata)

MicroRiPcor <- ggplot(studyFinal$Metadata, aes(AdjMicroglia, AdjustedRiP, color = Cohort))+
  theme_classic() +
  labs(x = "Microglia MGP (Adjusted)", y = "log(Normalized RiP) (Adjusted)") +
  geom_point(show.legend = F) +
  scale_color_manual(values = MoviePalettes$MadMaxDesert[c(1,7)]) +
  annotate("text", x = 0.3, y = 16.65, label = paste0("Cor =  ",
                                                     signif(MicroRiPCorStat$estimate, digits = 2),
                                                     ", p = ",
                                                     signif(MicroRiPCorStat$p.value, digits = 2)))

ggarrange(ggarrange(AdjMicrogliaPlot, AdjRiPPlot, nrow = 2), MicroRiPcor, nrow = 1)
ggsave(paste0(ResultsPath, "MicrogliaRiPcor.pdf"), device = "pdf", width = 9, height = 5, useDingbats = F)


ChIPResults %>% filter(symbol %in% (DESeqResults3 %>% filter(padj < 0.05, log2FoldChange > 0) %>% .$gene_name)) %>%
  group_by(region_type) %>% summarise(n()) %>% data.frame()

ChIPResults %>% filter(symbol %in% (DESeqResults3 %>% filter(padj < 0.05, log2FoldChange < 0) %>% .$gene_name)) %>%
  group_by(region_type) %>% summarise(n()) %>% data.frame()


GeneLength <- sapply((annoFileCollapsed %>% data.frame %>% filter(!is.na(symbol),
                                                                  !gene_type %in%  c("snoRNA", "misc_RNA",
                                                                                     "snRNA", "miRNA")) %>%
                                                                    .$symbol %>% unique), function(gene){
  Data = annoFileCollapsed %>% data.frame %>% filter(symbol == gene, region_type != "Up1to5Kb") %>% arrange(start)
  data.frame(Symbol = gene, Length = Data$end[nrow(Data)] - Data$start[1])
}, simplify = F) %>% rbindlist %>% data.frame

GeneLength$GeneType <- annoFileCollapsed$gene_type[match(GeneLength$Symbol, annoFileCollapsed$symbol)]

GeneLength %<>% filter(!GeneType %in%  c("snoRNA", "misc_RNA",
                                        "snRNA", "miRNA", "sense_intronic"))
GeneLength$DEgene <- "NS"
GeneLength$DEgene[GeneLength$Symbol %in% (DESeqResults3 %>% filter(padj < 0.05, log2FoldChange < 0) %>% .$gene_name)] <- "AgeDown"
GeneLength$DEgene[GeneLength$Symbol %in% (DESeqResults3 %>% filter(padj < 0.05, log2FoldChange > 0) %>% .$gene_name)] <- "AgeUp"

ggplot(GeneLength %>%  filter(Length > 5000), aes(log(Length))) + geom_density(aes(fill = DEgene, color = DEgene), alpha = 0.3)
