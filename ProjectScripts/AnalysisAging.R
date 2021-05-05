PeakLoc = "data/Peaks/" #This is the location (local) of the NarrowPeak files
CountMatrixLoc = "data/all_counts.tsv.gz" #This is the count matrix of our samples in our peaks
CellTypePeakCountLoc = "data/NeuN_peak_counts.tsv" #This is the count matrix of our samples is the cell type peaks

if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

################## Metadata ##############################################

source("ProjectScripts/PrepareMeta.R")

#######################Create files to analyse the called peaks ##############################

#Filter the samples 
Metadata %<>% filter(Age > 15) #this part is to filter samples
Metadata$Agef <- cut(Metadata$Age, breaks = 5, ordered_result = T)

AgeDist <- Metadata %>% group_by(Cohort2, Sex) %>%
  summarise(n = n(),
            Mean = mean(Age), Min = min(Age), Max = max(Age)) %>%
  data.frame()

AgeDist$Age = apply(AgeDist[,4:6], 1, function(x){
  paste0(round(x[1], digits = 0),
         " (", round(x[2], digits = 0),
         "-", round(x[3], digits =  0), ")")
})

write.table(AgeDist %>% select(Cohort2, Sex, n, Age), file = paste0(ResultsPath, "AgeDist.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

ggplot(Metadata, aes(FinalBatch, Age)) +
  theme_minimal() +
  labs(x = "") +
  geom_boxplot(outlier.shape = NA, aes(fill = Cohort2), alpha = 0.7) +
  geom_jitter(width = 0.2, aes(color = Sex)) +
  scale_color_manual(values = MoviePalettes$BugsLife[c(4,2)]) +
  scale_fill_manual(values = MoviePalettes$MoonRiseKingdomColors [c(8,2)], name = "Cohort")

ggsave(paste0(ResultsPath, "AgeCohortBatchDist.pdf"), device = "pdf",
       width = 6, height = 4, dpi = 300, useDingbats = F )
ggsave(paste0(ResultsPath, "AgeCohortBatchDist.jpg"), device = "jpg",
       width = 4, height = 2.5, dpi = 300)
closeDev()

# First for peaks called individually
InputPeakSingleCalled <- list()
InputPeakSingleCalled$PeakData <- as.list(as.character(Metadata$SampleID))
names(InputPeakSingleCalled$PeakData) <- make.names(as.character(Metadata$SampleID))

InputPeakSingleCalled$PeakData <- lapply(InputPeakSingleCalled$PeakData, function(sbj){
  fileName = paste0(PeakLoc, Metadata %>% filter(SampleID == sbj) %>% .$PeakFileSample)
  temp <- read.table(fileName, header = F, sep = "\t")
  names(temp) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")
  temp <- temp[!grepl("GL|hs", temp$CHR),] %>% droplevels()
  temp %<>% mutate(Length = END-START,
                   pValue = 10^(-pPvalue),
                   CHR = paste0("chr", CHR)) %>% filter(pValue < 10^(-7))
  temp %<>% arrange(CHR, START) %>% as(.,"GRanges")
  
  # #This is because the peaks were called based on hg19
  # temp2 <- LiftOver_hg19to38(temp)
  # temp2 %<>% data.frame
  # temp2$seqnames <- sapply(temp2$seqnames, function(x){
  #   gsub("chr", "", x)
  # })
  # 
  # temp2 %>% as(., "GRanges")
  
  temp
})#, mc.cores = detectCores())

InputPeakSingleCalled$UniquePeaks <- sapply(names(InputPeakSingleCalled$PeakData), function(sbj){
  temp <- InputPeakSingleCalled$PeakData[[sbj]]
  temp2 <- InputPeakSingleCalled$PeakData[!grepl(sbj, names(InputPeakSingleCalled$PeakData))]
  tempGRangesList <- GRangesList(temp2)
  OverlapPeaks <- findOverlaps(temp, tempGRangesList)
  temp[-queryHits(OverlapPeaks)]
}, simplify = FALSE)

InputPeakSingleCalled$SampleStat <- data.frame(SampleName = names(InputPeakSingleCalled$PeakData),
                                               activemotif_id = sapply(names(InputPeakSingleCalled$PeakData), function(x){
                                                 gsub("X", "", x)
                                               }),
                                               TotalPeaks = lapply(InputPeakSingleCalled$PeakData, function(sbj){
                                                 sbj %>% data.frame %>% nrow
                                               }) %>% unlist,
                                               TotalCoverage = lapply(InputPeakSingleCalled$PeakData, function(sbj){
                                                 sbj %>% data.frame() %>% .$width %>% sum
                                               }) %>% unlist,
                                               UniquePeaks = lapply(InputPeakSingleCalled$UniquePeaks, function(sbj){
                                                 sbj %>% data.frame %>% nrow
                                               }) %>% unlist,
                                               MedianPeakLengthAll = lapply(InputPeakSingleCalled$PeakData, function(sbj){
                                                 sbj %>% data.frame %>% .$width %>% median()
                                               }) %>% unlist,
                                               MedianPeakLengthUnique = lapply(InputPeakSingleCalled$UniquePeaks, function(sbj){
                                                 sbj %>% data.frame %>% .$width %>% median()
                                               }) %>% unlist)



InputPeakSingleCalled$SampleStat %<>% mutate(UniquePeakPercent = signif(100*UniquePeaks/TotalPeaks, digits = 3),
                                 TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))


InputPeakSingleCalled$SampleStat <- merge(InputPeakSingleCalled$SampleStat,  Metadata %>% select(-matches("File|seque|input|ip_|sample_id|comm", ignore.case = T)), by.x = "SampleName", by.y = "SampleID")

SampleStatMelted <- gather(InputPeakSingleCalled$SampleStat, key = "Measure_type", value = "Value", numReads, TotalPeaks, UniquePeakPercent, TotalCovPercent, MedianPeakLengthAll, MedianPeakLengthUnique)
SampleStatMelted$Measure_type <- factor(SampleStatMelted$Measure_type, levels = c("numReads", "TotalPeaks",  "TotalCovPercent", "MedianPeakLengthAll",  "MedianPeakLengthUnique", "UniquePeakPercent"))

Plot <- ggplot(SampleStatMelted, aes(Age, Value, color = FinalBatch)) +
  theme_classic(base_size = 14) +
  labs(x="", y="") +
  geom_point() +
  geom_smooth(method = "auto", color = "black") +
  facet_wrap(~Measure_type, scales = "free", ncol = 3)

ggplot(SampleStatMelted %>% filter(FinalBatch != "batch1"), aes(Age, Value, color = FinalBatch)) +
  theme_classic(base_size = 14) +
  labs(x="", y="") +
  geom_point() +
  geom_smooth(method = "auto", color = "black") +
  facet_wrap(~Measure_type + condition , scales = "free", ncol = 2)

ggsave(paste0(ResultsPath,"SingleCalledStat", Cohort, ".pdf"), plot = Plot,
       device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F)

lm(numReads ~ Agef + Sex + FinalBatch, data = InputPeakSingleCalled$SampleStat) %>% summary
lm(UniquePeakPercent ~ Agef + Sex + FinalBatch + numReads, data = InputPeakSingleCalled$SampleStat) %>% summary
lm(TotalCovPercent ~ Agef + Sex + FinalBatch + numReads , data = InputPeakSingleCalled$SampleStat) %>% summary


SubData <- InputPeakSingleCalled$SampleStat 
ModUniquePeaks <- lm(UniquePeakPercent ~ Agef + Sex + FinalBatch, data = SubData)
ModTotCov <- lm(TotalCovPercent ~ Agef + Sex + FinalBatch, data = SubData)
ModnmReads <- lm(numReads ~ Agef + Sex + FinalBatch, data = SubData)


SubData$AdjUniquePeakPercent <- ModelAdj(ModUniquePeaks, adj=data.frame(effect = c("FinalBatch", "Sex"), adjValue=c(0, 0)))
SubData$AdjTotCovPercent <- ModelAdj(ModTotCov, adj=data.frame(effect = c("FinalBatch", "Sex"), adjValue=c(0, 0)))
SubData$AdjnumReads <- ModelAdj(ModnmReads, adj=data.frame(effect = c("FinalBatch", "Sex"), adjValue=c(0, 0)))

SubDataMelt <- gather(SubData, key = "Measure_type", value = "Value", matches("Adj"))

Plot <- ggplot(SubDataMelt, aes(Age, log(Value), color = FinalBatch)) +
  theme_classic(base_size = 14) +
  theme_minimal() +
  labs(x="", y="log(value)") +
  geom_smooth(method = "auto", color = "black", size = 0.5) +
  scale_color_manual(values =  MoviePalettes$MoonRiseKingdomColors[c(2, 4, 6, 10)]) +
  geom_point(size = 1) +
  facet_wrap(~Measure_type, scales = "free", nrow = 1)
ggsave(paste0(ResultsPath,"SingleCalledStatBatchCorrected", Cohort, ".pdf"), plot = Plot, device = "pdf",
       width = 8, height = 2, dpi = 300, useDingbats = F )
closeDev()

################ Repeat for peaks called on all samples combined ############
InputPeakAllCalled <- list()
InputPeakAllCalled$PeakData <- read.table(paste0(PeakLoc, "all.narrowPeak.gz"), header = F, sep = "\t")
names(InputPeakAllCalled$PeakData) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")[1:ncol(InputPeakAllCalled$PeakData)]
InputPeakAllCalled$PeakData %<>% mutate(pValue = 10^(-pPvalue)) %>% filter(pValue < 10^(-7))
InputPeakAllCalled$PeakData <- InputPeakAllCalled$PeakData[!grepl("GL|hs", InputPeakAllCalled$PeakData$CHR),] %>% droplevels()
InputPeakAllCalled$PeakData %<>% arrange(CHR, START) %>% mutate(CHR = paste0("chr", CHR)) %>% as(.,"GRanges")

# #Lift over hg19 to hg38
# InputPeakAllCalled$PeakData <- LiftOver_hg19to38(InputPeakAllCalled$PeakData) %>% data.frame()
# InputPeakAllCalled$PeakData$seqnames <- sapply(InputPeakAllCalled$PeakData$seqnames, function(x){
#   gsub("chr", "", x)
# })
# 
# 
# InputPeakAllCalled$PeakData %<>% as(., "GRanges")

InputPeakAllCalled$Summary <- data.frame(TotalPeaks = InputPeakAllCalled$PeakData %>% data.frame %>% nrow,
                                         TotalCoverage = InputPeakAllCalled$PeakData %>% data.frame %>% .$width %>% sum)

InputPeakAllCalled$Summary %<>% mutate(TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))
InputPeakAllCalled$PeakData %>% data.frame() %>% .$width %>% summary()                                         
                                           

############################# HTseq counts ######################################################################
HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")

HTseqCounts %<>% select(-"X.data2.parkome.chipseq.aging_data_second_repeat.results.bamfiles.60_05H1_00ICHauk_H3K27Ac_hs_i88_ext_200.bam")
#names(HTseqCounts)[names(HTseqCounts) == "X.data2.parkome.chipseq.aging_data_second_repeat.results.bamfiles.60_05H1_00ICHauk_H3K27Ac_hs_i88_ext_200.bam"] <- "X60batch4"

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub(".data2.parkome.chipseq.aging_data_second_repeat.results.bamfiles.0?", "", x)
  strsplit(x,"_")[[1]][1]
})
HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

#Remove peaks in contig regions
HTseqCounts <- HTseqCounts[!grepl("GL|hs", HTseqCounts$CHR),] %>% droplevels()

#Blacklisted peaks
BlackListed <- read.table("data/H3K27Ac_black_list.bed", header = F, sep = "\t")
BlackListed %<>% mutate(Peak.Location = paste0("chr", V1, ":", V2, "-", "V3")) 
names(BlackListed)[1:3] <- c("CHR", "START", "END")

#Remove regions on Mitochondrial DNA
BlackListed <- BlackListed[!grepl("^M", BlackListed$CHR),] %>% droplevels()

#Filter the blacklisted peaks from the HTseq matrix
blackListedPeaks <- findOverlaps(query = HTseqCounts %>% as(., "GRanges"), subject = BlackListed %>% as(., "GRanges"), maxgap = 0, type = "any", select = "all")

HTseqCounts <- HTseqCounts[-queryHits(blackListedPeaks),]

#Remove peaks in contig regions and peaks with p > 10^-7 (as well as blacklisted peaks)
HTseqCounts <- HTseqCounts %>% filter(Geneid %in% as.character(InputPeakAllCalled$PeakData$PeakName)) %>% droplevels()


#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Length", as.character(Metadata$SampleID)))

# HTseqCounts <- merge(data.frame(InputPeakAllCalled$PeakData) %>%
#                        select(PeakName, seqnames, start, end, width),
#                      HTseqCounts, by.x = "PeakName", by.y = "Geneid", sort = F)
# 
# names(HTseqCounts)[1:5] <- c("Geneid", "CHR", "START", "END", "Length")

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("^X")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  

AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29", MetaSamleCol = "SampleID", countSampleRegEx = "^X",
                                     MetaCol = c("SampleID", "Sex", "Age", "Agef", "AgeGroup", "PMI", "Cohort", "FinalBatch", "numReads"))


##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
CellEstimateList <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")
CellEstimateList2 <- GetCellularProportions2(AllCalledData$SampleInfo,
                                             AllCalledData$countsMatrixAnnot,
                                             MetaSamplCol = "SampleID", normCol = "MeanRatioOrg")
CellEstimateList3 <- GetCellularProportions3(AllCalledData$SampleInfo,
                                             AllCalledData$countsMatrixAnnot,
                                             MetaSamplCol = "SampleID")


#Comparison to Human microglia expression from Galatro et al. 2017
MicrogliaAging <- read.table("data/MicrogliaHumanAgingGalatro.txt", header = T, sep = "\t")
MicogliaAll <- read.table("data/MicrogliaVsWholeBrainGalatro.txt", header = T, sep = "\t")

MicrogliaSpecificHuman <- MicogliaAll %>% filter(gliaVSbrain_logFC > 3, adj.P.Val < 0.05) %>% .$GeneSymbol
MicrogliaMSPgenes <- CellEstimateList$CellTypeSpecificPeaks$Microglia %>% filter(PeakName %in% rownames(CellEstimateList$PCAcellType$Microglia_MSP$rotation)) %>% .$symbol %>% unique
MicrogliaMSPgenes <- rownames(CellEstimateList3$PCAcellType$Microglia_Genes_MSP$rotation)

OverlapMicrogliaAll <- intersect(MicrogliaSpecificHuman, MicrogliaMSPgenes)
OverlapAgeUP <- intersect((MicrogliaAging %>% filter(logFC > 0) %>% .$GeneSymbol), MicrogliaMSPgenes)
OverlapAgeDown <- intersect((MicrogliaAging %>% filter(logFC < 0) %>% .$GeneSymbol), MicrogliaMSPgenes)



BackGroundBed <- sapply(names(CellEstimateList$CellTypeSpecificPeaks), function(cellTypeName){
  cellType = CellEstimateList$CellTypeSpecificPeaks[[cellTypeName]]
  sapply(unique(cellType$PeakName), function(Peak){
    data = cellType %>% filter(PeakName == Peak)
    
    if(nrow(data) > 1){
      data$symbol = paste0(data$symbol, collapse = "/")
      data = data[1,]
    }
    Background = c(500, 5*10^5)
    data.frame(Chr = data$BroadPeak.seqnames, Start = c(data$BroadPeak.start-Background, rep(data$BroadPeak.end, 2)),
               End = c(rep(data$BroadPeak.start, 2), data$BroadPeak.end + Background),
               PeakName = paste0(Peak, c("NarrowUp", "WideUP", "NarrowDown", "WideDown")),
               OrgPeak = Peak,
               GeneSymbol = data$symbol,
               cellType = cellTypeName) %>%
      mutate(Width = End - Start,
             UniqueCol = paste(PeakName, cellType,  sep = "_"))
  }, simplify = F) %>% rbindlist() %>% data.frame()
}, simplify = F) %>% rbindlist() %>% data.frame()

write.table(BackGroundBed, paste0(ResultsPath, "BackgroundCellPeaks.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

AllCalledData$SampleInfo <- CellEstimateList$Metadata

countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",CorMethod = "pearson",countSampleRegEx = "^X",MetaSamleCol = "SampleID", MetaSamleIDCol = "SampleID",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))

#Get the p-values for association of each covariates with the first 5 PCs
PCAsamples <- prcomp(t(countMatrixFullAllCalled$CPMdata), scale. = T)
countMatrixFullAllCalled$Metadata %<>% mutate(PC1 = PCAsamples$x[,1],
                                              PC2 = PCAsamples$x[,2],
                                              PC3 = PCAsamples$x[,3],
                                              PC4 = PCAsamples$x[,4],
                                              PC5 = PCAsamples$x[,5]) 
VarExplained <- PCAsamples %>% summary() %>% .$importance %>%
  .[2, 1:sum(grepl("^PC", names(countMatrixFullAllCalled$Metadata)))]*100 


CovarPvalues <- sapply(grep("^PC", names(countMatrixFullAllCalled$Metadata), value = T), function(PC){
  temp <- lm(as.formula(paste0(PC, "~ Age + Sex + FinalBatch + Oligo_MSP + Microglia_MSP + Endothelial_MSP+ NeuNall_MSP")),
             data = countMatrixFullAllCalled$Metadata) %>% summary
  temp$coefficients[-1,4]
}, simplify = F) %>% do.call(cbind, .) %>% data.frame()

names(CovarPvalues) <- paste0(names(CovarPvalues), "(", round(VarExplained, digits = 1), "%)")
CovarPvalues %<>% mutate(Variable = factor(rownames(CovarPvalues), levels = rownames(CovarPvalues)))

levels(CovarPvalues$Variable) <- c("Age", "Sex", grep("FinalBatch", levels(CovarPvalues$Variable), value = T),
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


countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|^X"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("^X")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)

lmMod1 <- lm(log(RiP_NormMeanRatioOrg)~Agef + Sex + FinalBatch + Oligo_MSP + NeuNall_MSP, data = countMatrixFullAllCalled$Metadata)
lmMod2 <- lm(log(RiP_NormMeanRatioOrg)~Age + Sex + FinalBatch + NeuNall_MSP + Oligo_MSP, data = countMatrixFullAllCalled$Metadata)

countMatrixFullAllCalled$Metadata$AdjustedRiP <- ModelAdj(lmMod2,
                                                          adj=data.frame(effect = c("FinalBatch", "Sex",
                                                                                    "NeuNall_MSP", "Oligo_MSP"), adjValue=c(0, 0, 0.5, 0.5)))
ggplot(countMatrixFullAllCalled$Metadata, aes(Age, AdjustedRiP, color = FinalBatch)) +
  theme_minimal()+
  theme(legend.position = "bottom") +
  labs(y = "Adjusted log(RiP)") +
  geom_point() + geom_smooth(color = "black") +
  scale_color_manual(values =  MoviePalettes$MoonRiseKingdomColors[c(2, 4, 6, 10)], name = "")


AgeChanges <- sapply(names(countMatrixFullAllCalled$Metadata)[grepl("_MSP",
                                                                    names(countMatrixFullAllCalled$Metadata))],
                     function(celltype){
  lm(as.formula(paste0(celltype, "~Age + Sex + FinalBatch")),
     data = countMatrixFullAllCalled$Metadata) %>% summary %>% .$coef
}, simplify = F)

AgeChangesDF <- sapply(names(AgeChanges), function(CellType){
  data <- AgeChanges[[CellType]] %>% data.frame()
  names(data)[4] <- "pValue"
  data[2,c(1,4)] %>% data.frame() %>% mutate(CellType = CellType)
}, simplify = F) %>% rbindlist()

write.table(AgeChangesDF, paste0(ResultsPath, "CellTypeAge.tsv"), sep = "\t",
            row.names = F, col.names = T)

#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


Model = as.formula(~Agef + Sex + FinalBatch + Oligo_MSP + NeuNall_MSP)
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")
                                  

DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.L_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.L") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.Q_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.Q") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.C_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.C") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")


ggMA(resultObject =  DESegResultsAge.L_FullAll %>% filter(!duplicated(PeakName)),
     lims = c(-1.2,1.2),
     geneColName = PeakName, contours = T, colour_pvalue = T )
ggsave(paste0(ResultsPath,"DiscoveryMAplot.pdf"), device = "pdf",
       width = 10, height = 4, dpi = 300, useDingbats = F)


ADgenesSims <- read.table("data/ADgeneSims2020.txt", header = T, sep = "\t")
ADgenesMeta <- read.table("data/ADsnps.txt", header = T, sep = "\t")

ADsignif <- DESegResultsAge.L_FullAll %>% filter(!duplicated(Peak_Gene), padj < 0.05) %>% arrange(padj) %>%
  select(PeakName, symbol, log2FoldChange, padj) %>%
  filter(symbol %in% ADgenesSims$GeneSymbol) %>% .$symbol %>% table() %>% data.frame()
names(ADsignif) <- c("GeneSymbol", "Number")

ADall <- DESegResultsAge.L_FullAll %>% filter(!duplicated(Peak_Gene), !is.na(padj)) %>% arrange(padj) %>%
  select(PeakName, symbol, log2FoldChange, padj) %>%
  filter(symbol %in% ADgenesSims$GeneSymbol) %>% .$symbol %>% table() %>%  data.frame()

names(ADall) <- c("GeneSymbol", "Number")

ADgenePeaks <- merge(ADsignif, ADall, by = "GeneSymbol",
                     all.x = T, all.y = T, suffixes = c("_Signif", "_All"))
ADgenePeaks$Number_Signif[is.na(ADgenePeaks$Number_Signif)] <- 0
ADgenePeaks %<>% arrange(desc(Number_Signif/Number_All), Number_All)

ADgenePeaks$Change <- sapply(ADgenePeaks$GeneSymbol, function(gene){
  temp <- DESegResultsAge.L_FullAll %>%
    filter(!duplicated(Peak_Gene), padj < 0.05, symbol == gene) %>%
    select(PeakName, symbol, log2FoldChange, padj)
  if(nrow(temp) == 0){
    NA
  } else {
    if(sum(temp$log2FoldChange < 0) == nrow(temp)){
      "Hypoacetylated"
    } else if(sum(temp$log2FoldChange > 0) == nrow(temp)){
      "Hyperacetylated"
    } else {
      "Mixed"
    }
  }
})

ADgenePeaks %<>% mutate("DARs (Peaks)" = paste0(Number_Signif, " (", Number_All, ")"))
ADgenePeaks <- merge(ADgenePeaks %>% select(-matches("Number")),
                     ADgenesSims, by = "GeneSymbol", sort = F)

write.table(ADgenePeaks, paste0(ResultsPath, "ADgeneDARs.tsv"), sep = "\t",
            row.names = F, col.names = T) 

TotalPeaks <- DESegResultsAge.L_FullAll %>%
  filter(!duplicated(PeakName), !is.na(padj)) %>% nrow

SignifPeaks <- DESegResultsAge.L_FullAll %>%
  filter(!duplicated(PeakName), !is.na(padj), padj < 0.05) %>% nrow

TotalADPeaks <- DESegResultsAge.L_FullAll %>%
  filter(!is.na(padj), symbol %in% ADgenesSims$GeneSymbol) %>% filter(!duplicated(PeakName)) %>% nrow

TotalSignifADPeaks <- DESegResultsAge.L_FullAll %>%
  filter(padj < 0.05, symbol %in% ADgenesSims$GeneSymbol) %>% filter(!duplicated(PeakName)) %>% nrow

phyper(TotalSignifADPeaks-1, TotalADPeaks,
       TotalPeaks - TotalADPeaks, SignifPeaks, lower.tail = F)


TotalGenes <- DESegResultsAge.L_FullAll %>%
  filter(!is.na(symbol), !is.na(padj)) %>% .[!grepl("^MIR|^SNO", .$symbol),] %>%
  .$symbol %>% unique %>% length

SignifGenes <- DESegResultsAge.L_FullAll %>%
  filter(!is.na(symbol), padj < 0.05) %>% .[!grepl("^MIR|^SNO", .$symbol),] %>%
  .$symbol %>% unique %>% length

TotalADGenes <- DESegResultsAge.L_FullAll %>%
  filter(!is.na(padj), symbol %in% ADgenesSims$GeneSymbol) %>%
  .$symbol %>% unique %>% length

SignifADGenes <- DESegResultsAge.L_FullAll %>%
  filter(padj < 0.05, symbol %in% ADgenesSims$GeneSymbol) %>%
  .$symbol %>% unique %>% length

phyper(SignifADGenes-1, TotalADGenes,
       TotalGenes - TotalADGenes, SignifADGenes, lower.tail = F)


#Looking at the association with peak length
PeakDist <- DESegResultsAge.L_FullAll %>% filter(!duplicated(PeakName)) %>%
  select(PeakName, Peak.width, baseMean, stat, pvalue, padj)
PeakDist$Change <- sapply(PeakDist$stat, function(x){
  if(x < 0){
    "Hypo"
  } else {
    "Hyper"
  }
})

PeakDist$baseCount <- sapply(PeakDist$baseMean, function(x){
  if(x < 100){
    "low"
  } else {
    "high"
  }
})

ggplot(PeakDist, aes(log10(Peak.width), -log10(pvalue))) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept = -log10(0.05),color = "red", lty = "dashed") +
  facet_wrap(~baseCount + Change, scales = "free")


#Enrichment analysis using GREAT
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


AgingEnrich <- GetChIPenrich(DESegResultsAge.L_FullAll)

save.image(paste0(ResultsPath, "WS_", Cohort, ".Rda"))
saveRDS(DESegResultsAge.L_FullAll, paste0(ResultsPath, "DESegResultsAge.L_FullAll.Rds"))
saveRDS(DESeqOutAll_Full, paste0(ResultsPath, "DESeqOutAll_Full.Rds"))

