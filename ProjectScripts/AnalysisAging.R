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
AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")


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
lmMod2 <- lm(log(RiP_NormMeanRatioOrg)~Agef + Sex + FinalBatch + NeuNall_MSP + Oligo_MSP + Microglia_MSP, data = countMatrixFullAllCalled$Metadata)
lmMod3 <- lm(log(RiP_NormMeanRatioOrg)~Agef + Sex + FinalBatch + NeuNall_MSP + Oligo_MSP + Endothelial_MSP + Microglia_MSP, data = countMatrixFullAllCalled$Metadata)

#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


Model = as.formula(~Agef + Sex + FinalBatch + Oligo_MSP + NeuNall_MSP + Microglia_MSP)
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")
                                  

DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.L_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.L") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.Q_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.Q") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.C_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.C") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")



save.image(paste0(ResultsPath, "WS_", Cohort, ".Rda"))
saveRDS(DESegResultsAge.L_FullAll, paste0(ResultsPath, "DESegResultsAge.L_FullAll.Rds"))
saveRDS(DESeqOutAll_Full, paste0(ResultsPath, "DESeqOutAll_Full.Rds"))

