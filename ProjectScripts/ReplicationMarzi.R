source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
packageF("tabulizer")

ResultsPath = "MarziAD_OurPeaks"

CountMatrixLoc = "data/all_counts.tsv.gz" #This is the count matrix of our samples in our peaks
CellTypePeakCountLoc = "data/NeuN_peak_counts.tsv" #This is the count matrix of our samples is the cell type peaks

CellTypePeakCountLoc = "DataMarzi/marzi_counts.tsv" #Count matrix of data from Marzi et al in the cell type peaks
#CountMatrixLoc = "DataMarzi/MarziPaperCounts.tsv"  #Count matrix from Marzi et al 
CountMatrixLoc = "DataMarzi/marzi_reads_on_our_aging_peaks.tsv.gz"  #Count matrix of data from Marzi et all quantified in our peaks
Cohort = "MarziAD"

if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")

ResultsDiscovery <- readRDS("AgingResults/DESegResultsAge.L_FullAll.Rds")
Deseq2OutDiscovery <- readRDS("AgingResults/DESeqOutAll_Full.Rds")

################## Metadata ##############################################
MetadataSup <- extract_tables("DataMarzi//41593_2018_253_MOESM1_ESM.pdf", pages = c(27:30)) %>% lapply(data.frame) %>% rbindlist() #Reading the metadata from the supplementary table 
names(MetadataSup) <- c("SampleID", "Group", "BraakStage", "Age", "Sex", "NeuralProportion", "PMI", "Experiment")
MetadataSup <- MetadataSup[-c(1:3),]
MetadataSup$Age <- as.numeric(as.character(MetadataSup$Age))
MetadataSup$PMI <- round(as.numeric(as.character(MetadataSup$PMI))/60, digits = 1)
MetadataSup$NeuralProportionNumeric <- sapply(MetadataSup$NeuralProportion, function(x){
  if(grepl("%", x)){
    x = gsub("%", "", x) %>% as.character() %>% as.numeric()
    x/100
  } else {
    NA
  }
})

MetadataSup %<>% mutate(MergeColumn = paste(Group, Age, Sex, NeuralProportionNumeric, sep = "_"))

#Getting additional metadata
softDown("GSE102538", file = "Data/GSE102538.soft")
Metadata <- ReadSoft("Data/GSE102538.soft") %>% select(-matches("antib|Sample_source|Sample_platform|orga")) %>% data.frame()
names(Metadata) <- c("SampleID", "SampleName", "Group", "Age", "Sex", "CETS")
Metadata$Age <- as.numeric(as.character(Metadata$Age))
Metadata$CETS <- as.numeric(as.character(Metadata$CETS))

Metadata %<>% arrange(CETS)
Metadata$Eno2 <- c(0.65, 0.99, 0.56, 0.48, 0.56, 0.55, 0.69, 0.26, 1.27, 0.88,
                   0.42, 0.51, 2.07, 0.84, 0.52, 1.30, 0.90, 0.80, 0.19, 1.16,
                   0.34, 0.69, 0.34, 0.77, 1.43, 0.72, 0.33, 0.36, 0.29, 0.59,
                   0.42, 0.83, 0.49, 1.53, 0.98, 1.31, 0.85, 1.61, 0.65, 1.11,
                   0.72, 1.26, 0.80, 0.48,  2.28, 1.17, NA)
Metadata$deltaCTEno2 <- sapply(Metadata$Eno2, function(x){
  if(!is.na(x)){
    -log2(x) 
  } else {
    NA
  }
}) 

Metadata$Group <- sapply(as.character(Metadata$Group), function(x) gsub("C", "Control", x))
Metadata$Group <- factor(Metadata$Group, levels = c("Control", "AD"))

Metadata$NeuralProportion <- round(Metadata$CETS, digits = 2)
Metadata %<>% mutate(MergeColumn = paste(Group, Age, Sex, NeuralProportion, sep = "_"))

Metadata <- merge(Metadata %>% select(-NeuralProportion), MetadataSup %>% select(BraakStage, PMI, MergeColumn), by = "MergeColumn", all.x = T, sort = F )

Metadata %<>% select(-MergeColumn)
Metadata %<>% filter(!duplicated(SampleID))

Metadata$Agef <- cut(Metadata$Age, breaks = 3, ordered_result = T)

      
rownames(Metadata) <- Metadata$SampleID %>% as.character()
Metadata$FinalBatch <- NA
Metadata$Cohort <- "Marzi"
############################# HTseq counts ######################################################################
HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")
HTseqCounts$Chr <- sapply(as.character(HTseqCounts$Chr), function(x){
  gsub("chr", "", x)
})

HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
#HTseqCounts %<>% mutate(Peak.Location = paste0("chr", CHR, ":", START, "-", END))

names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub("results.bamfiles.|_extended.bam", "", x)
  strsplit(x,"_")[[1]][1]
})

#Remove peaks in contig regions
HTseqCounts <- HTseqCounts[!grepl("GL|hs", HTseqCounts$CHR),] %>% droplevels()

#Blacklisted peaks
BlackListed <- read.table("data/H3K27Ac_black_list.bed", header = F, sep = "\t")
BlackListed %<>% mutate(Peak.Location = paste0("chr", V1, ":", V2, "-", V3)) 
names(BlackListed)[1:3] <- c("CHR", "START", "END")

#Remove regions on Mitochondrial DNA
BlackListed <- BlackListed[!grepl("^M", BlackListed$CHR),] %>% droplevels()

#Filter the blacklisted peaks from the HTseq matrix
blackListedPeaks <- findOverlaps(query = HTseqCounts %>% as(., "GRanges"), subject = BlackListed %>% as(., "GRanges"), maxgap = 0, type = "any", select = "all")

HTseqCounts <- HTseqCounts[-queryHits(blackListedPeaks),]


#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", as.character(Metadata$SampleID)))

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("GSM")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  

#Keep only the peaks included in the discovery analysis (the peaks were filtered based on MACS2 pvalue)
HTseqCounts %<>% filter(Geneid %in% ResultsDiscovery$PeakName)


AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29", MetaSamleCol = "SampleID",
                                     countSampleRegEx = "GSM", MetaCol = c("SampleID", "Sex", "Age", "Agef", "PMI", "BraakStage",  "PMI", "deltaCTEno2", "FinalBatch", "Group", "Cohort"))


##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")


countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)),
                                               collapseBy = "PeakName",CorMethod = "pearson",countSampleRegEx = "GSM",MetaSamleCol = "SampleID", MetaSamleIDCol = "SampleID",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))

lm(RiP_NormMeanRatioOrg~Age+Group+Sex+Oligo_MSP+Endothelial_MSP, data = countMatrixFullAllCalled$Metadata) %>% summary()

#Filter peaks with low counts

countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|GSM"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("GSM")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("GSM")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("GSM")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)

#Get the pvalues for associasion of each covariate with the first 3 PCs
PCAsamples <- prcomp(t(countMatrixFullAllCalled$CPMdata[rowSums(countMatrixFullAllCalled$CPMdata) > 0,]), scale. = T)
countMatrixFullAllCalled$Metadata %<>% mutate(PC1 = PCAsamples$x[,1],
                                              PC2 = PCAsamples$x[,2],
                                              PC3 = PCAsamples$x[,3],
                                              PC4 = PCAsamples$x[,4],
                                              PC5 = PCAsamples$x[,5]) 
VarExplained <- PCAsamples %>% summary() %>% .$importance %>%
  .[2, 1:sum(grepl("^PC", names(countMatrixFullAllCalled$Metadata)))]*100 


CovarPvalues <- sapply(grep("^PC", names(countMatrixFullAllCalled$Metadata), value = T), function(PC){
  temp <- lm(as.formula(paste0(PC, "~ Age + Sex + Group + Oligo_MSP + Microglia_MSP + Endothelial_MSP+ NeuNall_MSP")),
             data = countMatrixFullAllCalled$Metadata) %>% summary
  temp$coefficients[-1,4]
}, simplify = F) %>% do.call(cbind, .) %>% data.frame()

names(CovarPvalues) <- paste0(names(CovarPvalues), "(", round(VarExplained, digits = 1), "%)")
CovarPvalues %<>% mutate(Variable = factor(rownames(CovarPvalues), levels = rownames(CovarPvalues)))

levels(CovarPvalues$Variable) <- c("Age", "Sex", "Group",
                                   grep("MSP", levels(CovarPvalues$Variable), value = T))

CovarPvaluesMelt <- gather(CovarPvalues, key = "PC", value = "pValue", -Variable)
CovarPvaluesMelt %<>% mutate(pPvalue = -log10(pValue)) 


Plot  <- ggplot(CovarPvaluesMelt, aes(PC, Variable)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) +
  labs(x = "", y = "", title = "Association of variables with main PCs") +
  geom_tile(aes(fill = pPvalue), colour = "white") +
  scale_fill_gradient(high = "steelblue", low = "gray94", name = "-log10(p)") +
  geom_text(aes(label = signif(pValue, 2))) 
ggsave(paste0("AssociationWithPCs", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()



#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


Model = as.formula(~Agef + Sex + Group + Oligo_MSP + NeuNall_MSP + Microglia_MSP)
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")


DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
#DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

temp <- merge(ResultsDiscovery %>% filter(!duplicated(PeakName)) %>% select(PeakName, log2FoldChange, stat, pvalue, padj),
              DESegResultsGroup_FullAll %>% filter(!duplicated(PeakName)) %>% select(PeakName, log2FoldChange, stat, pvalue, padj),
              by = "PeakName", suffixes = c("_AgingAll", "_MarziAD"))

ggplot(temp, aes(stat_AgingAll, stat_MarziAD)) + geom_point()

cor.test(~stat_AgingAll + stat_MarziAD, data = temp %>% filter(pvalue_AgingAll < 0.05))

temp %>% filter(padj_AgingAll < 0.05, padj_MarziAD < 0.05) %>% dim()
temp %>% filter(padj_AgingAll < 0.05) %>% dim()
temp %>% filter(padj_MarziAD < 0.05) %>% dim()

dhyper(90, 1840, 615890-1840, 2689)


############################################################
#######  Repeat for aging, control samples only
############################################################
ResultsPath = "MarziAD_OurPeaks"
Cohort = "MarziAD"

if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

Metadata %<>% filter(Group == "Control")
Metadata$Agef2 <- cut(Metadata$Age, breaks = 5, ordered_result = T)
HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")
HTseqCounts$Chr <- sapply(as.character(HTseqCounts$Chr), function(x){
  gsub("chr", "", x)
})

HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
#HTseqCounts %<>% mutate(Peak.Location = paste0("chr", CHR, ":", START, "-", END))

names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub("results.bamfiles.|_extended.bam", "", x)
  strsplit(x,"_")[[1]][1]
})

#Remove peaks in contig regions
HTseqCounts <- HTseqCounts[!grepl("GL|hs", HTseqCounts$CHR),] %>% droplevels()

#Blacklisted peaks
BlackListed <- read.table("data/H3K27Ac_black_list.bed", header = F, sep = "\t")
BlackListed %<>% mutate(Peak.Location = paste0("chr", V1, ":", V2, "-", V3)) 
names(BlackListed)[1:3] <- c("CHR", "START", "END")

#Remove regions on Mitochondrial DNA
BlackListed <- BlackListed[!grepl("^M", BlackListed$CHR),] %>% droplevels()

#Filter the blacklisted peaks from the HTseq matrix
blackListedPeaks <- findOverlaps(query = HTseqCounts %>% as(., "GRanges"), subject = BlackListed %>% as(., "GRanges"), maxgap = 0, type = "any", select = "all")

HTseqCounts <- HTseqCounts[-queryHits(blackListedPeaks),]


#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", as.character(Metadata$SampleID)))

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("GSM")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  

#Keep only the peaks included in the discovery analysis (the peaks were filtered based on MACS2 pvalue)
HTseqCounts %<>% filter(Geneid %in% ResultsDiscovery$PeakName)


AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29", MetaSamleCol = "SampleID",
                                     countSampleRegEx = "GSM", MetaCol = c("SampleID", "Sex", "Age", "Agef", "Agef2", "PMI", "BraakStage",  "PMI", "deltaCTEno2", "FinalBatch", "Group", "Cohort"))


##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")


countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)),
                                               collapseBy = "PeakName",CorMethod = "pearson",countSampleRegEx = "GSM",MetaSamleCol = "SampleID", MetaSamleIDCol = "SampleID",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))

lm(RiP_NormMeanRatioOrg~Age+Sex+Oligo_MSP+Endothelial_MSP, data = countMatrixFullAllCalled$Metadata) %>% summary()

#Filter peaks with low counts

countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|GSM"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("GSM")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("GSM")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("GSM")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)

#Get the pvalues for associasion of each covariate with the first 3 PCs
PCAsamples <- prcomp(t(countMatrixFullAllCalled$CPMdata[rowSums(countMatrixFullAllCalled$CPMdata) > 0,]), scale. = T)
countMatrixFullAllCalled$Metadata %<>% mutate(PC1 = PCAsamples$x[,1],
                                              PC2 = PCAsamples$x[,2],
                                              PC3 = PCAsamples$x[,3],
                                              PC4 = PCAsamples$x[,4],
                                              PC5 = PCAsamples$x[,5]) 
VarExplained <- PCAsamples %>% summary() %>% .$importance %>%
  .[2, 1:sum(grepl("^PC", names(countMatrixFullAllCalled$Metadata)))]*100 


CovarPvalues <- sapply(grep("^PC", names(countMatrixFullAllCalled$Metadata), value = T), function(PC){
  temp <- lm(as.formula(paste0(PC, "~ Age + Sex + Oligo_MSP + Microglia_MSP + Endothelial_MSP+ NeuNall_MSP")),
             data = countMatrixFullAllCalled$Metadata) %>% summary
  temp$coefficients[-1,4]
}, simplify = F) %>% do.call(cbind, .) %>% data.frame()

names(CovarPvalues) <- paste0(names(CovarPvalues), "(", round(VarExplained, digits = 1), "%)")
CovarPvalues %<>% mutate(Variable = factor(rownames(CovarPvalues), levels = rownames(CovarPvalues)))

levels(CovarPvalues$Variable) <- c("Age", "Sex", "Group",
                                   grep("MSP", levels(CovarPvalues$Variable), value = T))

CovarPvaluesMelt <- gather(CovarPvalues, key = "PC", value = "pValue", -Variable)
CovarPvaluesMelt %<>% mutate(pPvalue = -log10(pValue)) 


Plot  <- ggplot(CovarPvaluesMelt, aes(PC, Variable)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) +
  labs(x = "", y = "", title = "Association of variables with main PCs") +
  geom_tile(aes(fill = pPvalue), colour = "white") +
  scale_fill_gradient(high = "steelblue", low = "gray94", name = "-log10(p)") +
  geom_text(aes(label = signif(pValue, 2))) 
ggsave(paste0("AssociationWithPCs", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()



#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


Model = as.formula(~Agef2 + Sex + Oligo_MSP + Endothelial_MSP)
DESeqOutAll_Full2 <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")


DESegResultsSex_FullAll2 <- GetDESeqResults(DESeqOutAll_Full2, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
#DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

DESegResultsAge.L_FullAll <- GetDESeqResults(DESeqOutAll_Full2, coef = "Agef2.L") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.Q_FullAll <- GetDESeqResults(DESeqOutAll_Full2, coef = "Agef2.Q") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.C_FullAll <- GetDESeqResults(DESeqOutAll_Full2, coef = "Agef2.C") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")


DicovRepResult <- merge(ResultsDiscovery %>% filter(!duplicated(PeakName)) %>% select(PeakName, stat, pvalue, padj, log2FoldChange),
                        DESegResultsAge.L_FullAll %>% filter(!duplicated(PeakName)) %>% select(PeakName, stat, pvalue, padj, log2FoldChange),
                        by = "PeakName", suffixes = c("_Discov", "_Replic"))

AllThreeResults <- merge(DicovRepResult, DESegResultsGroup_FullAll %>% filter(!duplicated(PeakName)) %>% select(PeakName, stat, pvalue, padj, log2FoldChange),
                        by = "PeakName")

ggplot(AllThreeResults %>% filter( pvalue_Discov < 0.05), aes(stat_Discov, stat_Replic)) + geom_point()
cor.test(~stat_Discov + stat, data = AllThreeResults  %>% filter(pvalue_Discov < 0.05))

ggplot(AllThreeResults, aes(stat_Replic, stat)) + geom_point()


packageF("metap")
MetaP <- apply(DicovRepResult %>% select(matches("pvalue")), 1, function(x){
  allmetap(x, method = "all") %>% select(p) %>% t
}) %>% rbindlist() %>% data.frame()
names(MetaP) <- rownames(allmetap(c(4.223084e-09, 7.586691e-04), method = "all"))



DicovRepResult <- cbind(DicovRepResult, MetaP)

DicovRepResult$sumlogAdj <- p.adjust(DicovRepResult$sumlog, method = "BH")


ParkOmeSignif <- read.table("SignifPeaksOligEndoCorrect.tsv", header = T, sep = "\t")
ParkOmeSignifPos <- as.character(ParkOmeSignif %>%
                                   filter(!is.na(ParkOmeSignif$symbol), log2FoldChange > 0) %>% .$symbol)
ParkOmeSignifNeg <- as.character(ParkOmeSignif %>%
                                   filter(!is.na(ParkOmeSignif$symbol), log2FoldChange < 0) %>% .$symbol)


packageF("pROC")
Ranks <- which((DESegResultsAge.L_FullAll %>%
                 arrange(pvalue) %>%
                 filter(!duplicated(symbol)) %>%
                  filter(log2FoldChange > 0) %>%
                          .$symbol) %in% ParkOmeSignifPos)

temp <- DESegResultsAge.L_FullAll %>%
  arrange(pvalue) %>%
  filter(!duplicated(symbol)) %>% select(log2FoldChange, pvalue, symbol)

temp$SignifDiscovery <- sapply(temp$symbol, function(x){
  if(x %in% c(ParkOmeSignifPos, ParkOmeSignifNeg)){
    "Yes"
  } else {
    "No"
  }
})

plot.roc(roc(temp$SignifDiscovery, -log10(temp$pvalue)))

Ranks2 <- which((DESegResultsAge.L_FullAll %>%
                   arrange(pvalue) %>%
                   filter(!duplicated(symbol)) %>%
                   filter(log2FoldChange < 0) %>%
                   .$symbol) %in% ParkOmeSignifNeg)
                  
                

Ranks3 <- which((DESegResultsAge.L_FullAll %>% arrange(pvalue) %>%
                  filter(!duplicated(symbol)) %>% .$symbol) %in% as.character(ParkOmeSignif %>%
                                                                                filter(!is.na(ParkOmeSignif$symbol)) %>% .$symbol))


ADgenes <- read.table("data/ADgenes.txt", header = T, sep = "\t")
PDgenes <- read.table("data/PDgenes.txt", header = T, sep = "\t")
ALSgenes <- read.table("data/ALSgenes.txt", header = T, sep = "\t")

DiseaseGenes <- list(AD = DESegResultsAge.L_FullAll %>% filter(symbol %in% ADgenes$Gene, !is.na(padj)) %>% 
                       filter(!duplicated(Peak_Gene)) %>% select(symbol, padj) %>% mutate(Disease = "AD"),
                     PD = DESegResultsAge.L_FullAll %>% filter(symbol %in% PDgenes$Gene, !is.na(padj)) %>% 
                       filter(!duplicated(Peak_Gene)) %>% select(symbol, padj) %>% mutate(Disease = "PD"),
                     ALS = DESegResultsAge.L_FullAll %>% filter(symbol %in% ALSgenes$Gene, !is.na(padj)) %>% 
                       filter(!duplicated(Peak_Gene)) %>% select(symbol, padj) %>% mutate(Disease = "ALS"),
                     Other = DESegResultsAge.L_FullAll %>% filter(!symbol %in% c(as.character(ADgenes$Gene),
                                                                                 as.character(PDgenes$Gene),
                                                                                 as.character(ALSgenes$Gene)), !is.na(padj)) %>% 
                       filter(!duplicated(Peak_Gene)) %>% select(symbol, padj) %>% mutate(Disease = "Other")) %>%
  rbindlist() %>% data.frame() %>% arrange(padj)

ggplot(DiseaseGenes, aes(Disease, -log10(padj))) +
  geom_violin() + geom_boxplot(outlier.shape = NA, width = 0.1)

#DESegResultsAgeMiddle_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "AgeGroupMiddle") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
#DESegResultsAgeOld_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "AgeGroupgoOld") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

#save.image(paste0(ResultsPath, "WS_", Cohort, ".Rda"))