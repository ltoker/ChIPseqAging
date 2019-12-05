source("/data/Rprojects/GeneralScripts/generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
ResultsPath = "AgingResults"

PeakLoc = "data/Peaks/"
CountMatrixLoc = "data/all_counts.tsv.gz"
CellTypePeakCountLoc = "data/NeuN_peak_counts.tsv"

Cohort = "Aging"

if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")



################## Metadata ##############################################

source("ProjectScripts/PrepareMeta.R")

#######################Create files to analyse the called peaks ##############################

# First for peaks called individually
InputPeakSingleCalled <- list()
InputPeakSingleCalled$PeakData <- as.list(as.character(Metadata$SampleID))
names(InputPeakSingleCalled$PeakData) <- make.names(as.character(Metadata$SampleID))

InputPeakSingleCalled$PeakData <- mclapply(InputPeakSingleCalled$PeakData, function(sbj){
  fileName = paste0(PeakLoc, Metadata %>% filter(SampleID == sbj) %>% .$PeakFileSample)
  temp <- read.table(fileName, header = F, sep = "\t")
  names(temp) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")
  temp <- temp[!grepl("GL|hs", temp$CHR),] %>% droplevels()
  temp %<>% mutate(Length = END-START, pValue = 10^(-pPvalue)) %>% filter(pValue < 10^(-7))
  temp %>% arrange(CHR, START) %>% as(.,"GRanges")
}, mc.cores = detectCores())

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

ggsave(paste0(ResultsPath,"SingleCalledStat", Cohort, ".png"), plot = Plot, device = "png", width = 10, height = 6, dpi = 300)
ggsave(paste0(ResultsPath,"SingleCalledStat", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300 )

lm(numReads ~ Agef + Sex + FinalBatch, data = InputPeakSingleCalled$SampleStat %>% filter(SampleName != "X72")) %>% summary
lm(UniquePeakPercent ~ Agef + Sex + FinalBatch + numReads, data = InputPeakSingleCalled$SampleStat %>% filter(SampleName != "X72")) %>% summary
lm(TotalCovPercent ~ Agef + Sex + FinalBatch + numReads , data = InputPeakSingleCalled$SampleStat %>% filter(SampleName != "X72")) %>% summary


SubData <- InputPeakSingleCalled$SampleStat %>% filter(SampleName != "X72")
ModUniquePeaks <- lm(UniquePeakPercent ~ Agef + Sex + FinalBatch, data = SubData)
ModTotCov <- lm(TotalCovPercent ~ Agef + Sex + FinalBatch, data = SubData)

SubData$AdjUniquePeakPercent <- ModelAdj(ModUniquePeaks, adj=data.frame(effect = c("FinalBatch", "Sex"), adjValue=c(0, 0)))
SubData$AdjTotCovPercent <- ModelAdj(ModTotCov, adj=data.frame(effect = c("FinalBatch", "Sex"), adjValue=c(0, 0)))

SubDataMelt <- gather(SubData, key = "Measure_type", value = "Value", matches("Adj"))

ggplot(SubDataMelt, aes(Age, Value, color = FinalBatch)) +
  theme_classic(base_size = 14) +
  labs(x="", y="") +
  geom_point() +
  geom_smooth(method = "auto", color = "black") +
  facet_wrap(~Measure_type, scales = "free", nrow = 2)

################ Repeat for peaks called on all samples combined ############
InputPeakAllCalled <- list()
InputPeakAllCalled$PeakData <- read.table(paste0(PeakLoc, "all.narrowPeak.gz"), header = F, sep = "\t")
names(InputPeakAllCalled$PeakData) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")[1:ncol(InputPeakAllCalled$PeakData)]
InputPeakAllCalled$PeakData %<>% mutate(pValue = 10^(-pPvalue)) %>% filter(pValue < 10^(-7))
InputPeakAllCalled$PeakData <- InputPeakAllCalled$PeakData[!grepl("GL|hs", InputPeakAllCalled$PeakData$CHR),] %>% droplevels()
InputPeakAllCalled$PeakData %<>% arrange(CHR, START) %>% as(.,"GRanges")



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

#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", as.character(Metadata$SampleID)))

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

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("^X")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  

AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29", MetaSamleCol = "SampleID", countSampleRegEx = "^X", MetaCol = c("SampleID", "Sex", "Age", "Agef", "AgeGroup", "PMI", "Cohort", "FinalBatch", "numReads"))


closeDev()


##### Get relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
AllCalledData$SampleInfo <- GetCellularProportions(AllCalledData$SampleInfo, MetaSamplCol = "SampleID")

pdf(paste0(ResultsPath, "SampleCorAllPeaks", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",CorMethod = "pearson",countSampleRegEx = "^X",MetaSamleCol = "SampleID", MetaSamleIDCol = "SampleID",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))
closeDev()



countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|^X"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("^X")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)


################# Check the best normalization method ###############################

# ## Calculate normalization factors using RLE
# DEStemp <- DESeqDataSetFromMatrix(countData = countMatrixDF %>% data.frame %>% select(matches("^X")) %>% as.matrix, colData = countMatrixFullAllCalled$Metadata, design = Model)
# DEStemp <-  estimateSizeFactors(DEStemp)
# DEStemp$RiP_RLE <- DEStemp$TotalCount/DEStemp$sizeFactor
# 
# MeltedDataAll <- DEStemp@colData %>% data.frame %>% gather(key = "MeasureType", value = "Value", TotalCount, library_size, RiP_NormAllCount, RiP_NormBackground, RiP_RLE, RiP_NormMeanRatioOrg, RiP_NormMeanRatioAll)
# MeltedDataAll$MeasureType <- factor(MeltedDataAll$MeasureType, levels = c("TotalCount", "library_size",  "RiP_NormAllCount", "RiP_NormBackground", "RiP_RLE", "RiP_NormMeanRatioOrg", "RiP_NormMeanRatioAll"))
# levels(MeltedDataAll$MeasureType) <- c("RiP", "LibrarySize",  "RiP/LibrarySize", "RiP/RoP", "RiP/sizeFactor", "RiP/MeanRatio", "RiP/MeanRatio2")
# 
# #Add individual vallues
# xLabFun <- function(x) signif(as.numeric(as.character(x)), digits = 2)
# 
# 
# MeltedDataAll %<>% mutate(MeasureType2 = MeasureType)
# 
# levels(MeltedDataAll$MeasureType2) <- sapply(levels(MeltedDataAll$MeasureType2), function(Type){
#   temp <- StatNormMethod[[Type]] %>% summary()
#   paste0(Type, " (r=", signif(temp$r.squared^0.5, digits = 2), ", p=",  signif(temp$coefficients[2,4], digits = 2), ")")
# })
# 
# 
# ggplot(MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm)), aes(H3K27gapdh_Norm, Value, color = condition)) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   labs(y = "Counts", x = "WB, H3K27gapdh_normalized (Final)", title = "All samples Parkome") +
#   geom_point(alpha = 0.9) +
#   scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") + 
#   #geom_smooth(method = "lm", color = "black") +
#   scale_y_continuous(labels = xLabFun) +
#   facet_wrap(~MeasureType2, scales = "free_y")
# ggsave(paste0("WB_ChipSeqCor_FinalnormalizedToReplicates", Cohort, ".png"))
# 
# 
# 
# ggplot(MeltedDataAll %>% filter(MeasureType == "RiP/MeanRatio", batch != "H"), aes(age, Value, fill = condition)) +
#   theme_bw(base_size = 14) +
#   theme(panel.grid = element_blank()) +
#   labs(x = "Age", y = "Normalized RiP", title = "All samples") +
#   geom_smooth(method = "lm", aes(fill = condition, color = condition), alpha = 0.3, size = 0.2) +
#   geom_point(size = 2, aes(color = condition)) +
#   scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
#   scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
#   scale_y_continuous(labels = xLabFun) +
#   facet_wrap(~condition, scales = "free_x")
# ggsave(paste0("AgeRipCorrelation", Cohort, ".png"))
# 
# 
# 
# #Detect outliers
# MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
# Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


Model = as.formula(~AgeGroup + Sex + FinalBatch + Oligo_MSP + Microglia_MSP)
DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, MetaSamleCol = "SampleID",SampleNameCol = "SampleID",
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg",
                             FullModel = Model, test = "Wald", FitType = "local")
                                  

DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.L_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.L") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.Q_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.Q") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge.C_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Agef.C") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

DESegResultsAgeYoung_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "AgeGroupYoung") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAgeOld_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "AgeGroupOld") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
