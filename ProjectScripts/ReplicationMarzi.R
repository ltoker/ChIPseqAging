library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

source("ProjectScripts/ProjectFunctions.R")

CellTypePeakCountLoc = "DataMarzi/marzi_counts.tsv" #Count matrix of data from Marzi et al in the cell type peaks

#CountMatrixLoc = "DataMarzi/MarziPaperCounts.tsv"  #Count matrix from Marzi et al 
CountMatrixLoc = "DataMarzi/marzi_reads_on_our_aging_peaks.tsv.gz"  #Count matrix of data from Marzi et all quantified in our peaks

if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

ResultsDiscovery <- readRDS(paste0(DiscoveryResults, "/DESegResultsAge.L_FullAll.Rds"))
Deseq2OutDiscovery <- readRDS(paste0(DiscoveryResults, "/DESeqOutAll_Full.Rds"))

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

Metadata$Agef <- cut(Metadata$Age, breaks = 5, ordered_result = T)

      
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

lm(RiP_NormMeanRatioOrg~Age+Group+Sex+Oligo_MSP+Endothelial_MSP+Microglia_MSP+NeuNall_MSP, data = countMatrixFullAllCalled$Metadata) %>% summary()
lm(NeuNall_MSP~Age+Group+Sex+Oligo_MSP+Microglia_MSP, data = countMatrixFullAllCalled$Metadata) %>% summary()
lm(Microglia_MSP~Age+Group+Sex+NeuNall_MSP+Oligo_MSP, data = countMatrixFullAllCalled$Metadata) %>% summary()
lm(Oligo_MSP~Age+Group+Sex+NeuNall_MSP+Microglia_MSP, data = countMatrixFullAllCalled$Metadata) %>% summary()

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


Model = as.formula(~Agef + Sex + Group + Oligo_MSP + NeuNall_MSP + Microglia_MSP + Endothelial_MSP)
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")


DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
#DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")



############################################################
#######  Repeat for aging, control samples only
############################################################
Cohort = "Marzi_Control"


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

lm(RiP_NormMeanRatioOrg~Age+Sex+Oligo_MSP+Endothelial_MSP + Microglia_MSP + NeuNall_MSP , data = countMatrixFullAllCalled$Metadata) %>% summary()

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

levels(CovarPvalues$Variable) <- c("Age", "Sex",
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


Model = as.formula(~Agef + Sex + Oligo_MSP + NeuNall_MSP + Microglia_MSP)
DESeqOutAll_Full2 <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")


DESegResultsSex_FullAll2 <- GetDESeqResults(DESeqOutAll_Full2, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
#DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

DESegResultsAge.L_FullAll <- GetDESeqResults(DESeqOutAll_Full2, coef = "Agef.L") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.Q_FullAll <- GetDESeqResults(DESeqOutAll_Full2, coef = "Agef.Q") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.C_FullAll <- GetDESeqResults(DESeqOutAll_Full2, coef = "Agef.C") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")



saveRDS(list(DESeqOutMarziAD = DESeqOutAll_Full,
             DESegResultsAD = DESegResultsGroup_FullAll,
             DESeqOutMarziAging = DESeqOutAll_Full2,
             DESegResultsAgaingL = DESegResultsAge.L_FullAll), file = paste0(DiscoveryResults,"/OutputMarzi.Rds"))

DicovRepResult <- merge(ResultsDiscovery %>% filter(!duplicated(PeakName)) %>% select(PeakName, stat, pvalue, padj, log2FoldChange),
                        DESegResultsAge.L_FullAll %>% filter(!duplicated(PeakName)) %>% select(PeakName, stat, pvalue, padj, log2FoldChange),
                        by = "PeakName", suffixes = c("_Discov", "_Replic"), all.x = T, all.y = T)

AllThreeResults <- merge(DicovRepResult, DESegResultsGroup_FullAll %>% filter(!duplicated(PeakName)) %>% select(PeakName, stat, pvalue, padj, log2FoldChange),
                        by = "PeakName", all.x = T, all.y = T)

ggplot(AllThreeResults %>% filter(pvalue_Discov < 0.05), aes(stat_Discov, stat_Replic)) + geom_point()
cor.test(~stat_Discov + stat_Replic, data = AllThreeResults  %>% filter(pvalue_Discov < 0.05))

GetHypergeometric <- function(DF = AllThreeResults, Col1, Col2) {
  Data = DF %>% filter(!(is.na(.data[[Col1]]) | is.na(.data[[Col2]])))
  AllRelPeaks = Data %>% nrow()
  SignifBoth = Data %>% filter(.data[[Col1]] < 0.05, .data[[Col2]] < 0.05) %>% nrow()
  Signif1 = Data %>% filter(.data[[Col1]] < 0.05) %>% nrow()
  Signif2 = Data %>% filter(.data[[Col2]] < 0.05) %>% nrow()
  
  phyper(SignifBoth-1, Signif1, AllRelPeaks-Signif1, Signif2, lower.tail = F)
}

GetHypergeometric(Col1 = "padj_Discov", Col2 = "padj")
GetHypergeometric(Col1 = "padj_Discov", Col2 = "padj_Replic")
GetHypergeometric(Col1 = "padj", Col2 = "padj_Replic")



           