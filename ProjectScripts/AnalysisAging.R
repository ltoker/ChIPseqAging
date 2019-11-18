source("/data/Rprojects/GeneralScripts/generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
ResultsPath = "AgingResults"

if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

#This files contains all the relevant file locations and filters samples relevant to the analysis. Comment the irrelevant lines
source("ProjectScripts/ConfigFile.R") 
source("ProjectScripts/PrepareMeta.R")

MetaChipSeq %<>% mutate(bamReads = paste0(dataDir, ChipSampleID, "_final.bam"),
                        bamControl = paste0(dataDir, ChipControlID, "_final.bam"),
                        Peaks = paste0(dataDir, "Peaks/", ChipSampleID, ".narrowPeakClean.gz"),
                        PeakCaller = "narrow")


plotMA = DESeq2::plotMA
Metadata %<>% mutate(SampleID = paste0("X", activemotif_id))

#Filter only the relevant samples
MetaChipSeq %<>% filter(SampleID %in% Metadata$activemotif_id) %>% droplevels()
annoFileCollapsed <- GetGenomeAnno(genome = "hg19")


################## Metadata ##############################################

MetaChipSeqGroupCalled <- MetaChipSeq
MetaChipSeqGroupCalled$Peaks <- as.character(MetaChipSeqGroupCalled$Peaks)
MetaChipSeqGroupCalled$Peaks[grepl("Case", MetaChipSeqGroupCalled$Condition)] <- CaseNarrowPeakLoc
MetaChipSeqGroupCalled$Peaks[grepl("Control", MetaChipSeqGroupCalled$Condition)] <- ControlNarrowPeakLoc

MetaChipSeqAllcalled <- MetaChipSeq

MetaChipSeqAllcalled$Peaks <- CombinedNarrowPeakLoc

###### #Get the peaks ########################################################################
SingleSampleCalled <- dba(sampleSheet = MetaChipSeq, minOverlap=1,
                          config=data.frame(RunParallel=TRUE, reportInit="SingleSampleCalled", DataType=DBA_DATA_GRANGES,
                                            AnalysisMethod="chip_DESEQ2", minQCth=30, fragmentSize=172,
                                            bCorPlot=TRUE, th=0.05, bUsePval=FALSE))

AllCalled <- dba(sampleSheet = MetaChipSeqAllcalled, minOverlap=1,
                 config=data.frame(RunParallel=TRUE, reportInit="GroupCalled", DataType=DBA_DATA_GRANGES,
                                   AnalysisMethod="chip_DESEQ2", minQCth=30, fragmentSize=200,
                                   bCorPlot=TRUE, th=0.05, bUsePval=FALSE))

#######################Create files to analyse the called peaks ##############################

# First for peaks called individually
InputPeakSingleCalled <- SingleSampleCalled[c("peaks", "merged", "called")]

names(InputPeakSingleCalled$peaks) <- paste0("X",  SingleSampleCalled$samples$SampleID)

InputPeakSingleCalled$merged %<>% data.frame %>% mutate(Length = END-START)

InputPeakSingleCalled$called <- cbind(InputPeakSingleCalled$merged, data.frame(InputPeakSingleCalled$called))
names(InputPeakSingleCalled$called)[grepl("^X", names(InputPeakSingleCalled$called))] <- paste0("X",  SingleSampleCalled$samples$SampleID)

InputPeakSingleCalled$called$SamplesCalled <- apply(InputPeakSingleCalled$called %>% select(matches("^X")), 1, sum)

InputPeakSingleCalled$called$PeakNum = rownames(InputPeakSingleCalled$called)

InputPeakSingleCalled$NarrowPeak <- as.list(names(InputPeakSingleCalled$called %>% select(matches("^X"))))
names(InputPeakSingleCalled$NarrowPeak) <- names(InputPeakSingleCalled$called %>% select(matches("^X")))

InputPeakSingleCalled$NarrowPeak <- mclapply(InputPeakSingleCalled$NarrowPeak, function(sbj){
  sbj = gsub("X", "", sbj)
  fileName = MetaChipSeq %>% filter(SampleID == sbj) %>% .$Peaks
  temp <- read.table(fileName, header = F, sep = "\t")
  names(temp) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")
  temp %<>% mutate(Length = END-START, pValue = 10^(-pPvalue))
  temp %>% arrange(CHR, START)
}, mc.cores = detectCores())


InputPeakSingleCalled$SampleStat <- data.frame(SampleName = names(InputPeakSingleCalled$NarrowPeak),
                                               activemotif_id = sapply(names(InputPeakSingleCalled$NarrowPeak), function(x){
                                                 gsub("X", "", x)
                                               }),
                                               TotalPeaks = lapply(InputPeakSingleCalled$peaks, function(sbj){
                                                 sbj %>% .[!grepl("GL|hs|MT", .$V1),] %>% nrow
                                               }) %>% unlist,
                                               TotalCoverage = lapply(InputPeakSingleCalled$NarrowPeak, function(sbj){
                                                 sbj %>% .[!grepl("GL|hs|MT", .$CHR),] %>% .$Length %>% sum
                                               }) %>% unlist,
                                               TotalPeaksMerged = InputPeakSingleCalled$called %>% .[!grepl("GL|hs|MT", .$CHR),] %>% select(matches("^X")) %>% apply(2, sum),
                                               UniquePeaks = InputPeakSingleCalled$called %>% filter(SamplesCalled == 1) %>% .[!grepl("GL|hs|MT", .$CHR),] %>% select(matches("^X")) %>% apply(2, sum) %>% unlist)



InputPeakSingleCalled$SampleStat %<>% mutate(UniquePeakPercent = signif(100*UniquePeaks/TotalPeaks, digits = 3),
                                 TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))


InputPeakSingleCalled$SampleStat <- merge(InputPeakSingleCalled$SampleStat,  Metadata %>% select(activemotif_id, rin, condition, sex, age, batch,library_size,
                                                                                                 Oligo_Genes, Astrocyte_Genes, Microglia_Genes, NeuronAll_Genes,
                                                                                                 h3k27_gapdh, h3k27_h3, H3K27gapdh, H3K27gapdh_Norm), by = "activemotif_id")

SampleStatMelted <- gather(InputPeakSingleCalled$SampleStat, key = "Measure_type", value = "Value", UniquePeakPercent, TotalCovPercent)
ggplot(SampleStatMelted, aes(age, Value)) +
  theme_bw(base_size = 14) +
  geom_point(aes(color = condition)) +
  scale_color_manual(values = c("cornflowerblue", "coral4", "orange", "chartreuse4", "darkgrey", "darkmagenta")) +
  facet_grid(Measure_type~condition, scales = "free_y")

ggsave(paste0("SingleCalledStat", Cohort, ".png"))


lm(UniquePeakPercent ~ condition + batch  + age, data = InputPeakSingleCalled$SampleStat) %>% summary
lm(TotalCovPercent ~ condition + batch + age, data = InputPeakSingleCalled$SampleStat) %>% summary


################ Repeat for peaks called on all samples combined ############
InputPeakAllCalled <- AllCalled[c("peaks", "merged", "called")]

names(InputPeakAllCalled$peaks) <- paste0("X",  AllCalled$samples$SampleID)

InputPeakAllCalled$merged %<>% data.frame %>% mutate(Length = END-START)

InputPeakAllCalled$called <- cbind(InputPeakAllCalled$merged, data.frame(InputPeakAllCalled$called))
names(InputPeakAllCalled$called)[grepl("^X", names(InputPeakAllCalled$called))] <- paste0("X", AllCalled$samples$SampleID)

InputPeakAllCalled$called$SamplesCalled <- apply(InputPeakAllCalled$called %>% select(matches("^X")), 1, sum)

InputPeakAllCalled$called$PeakNum = rownames(InputPeakAllCalled$called)

InputPeakAllCalled$NarrowPeak <- as.list(unique(MetaChipSeqAllcalled$Peaks))
names(InputPeakAllCalled$NarrowPeak) <- "All"

InputPeakAllCalled$NarrowPeak <- mclapply(InputPeakAllCalled$NarrowPeak, function(grp){
  temp <- read.table(grp, header = F, sep = "\t")
  names(temp) <- c("CHR", "START", "END", "PeakName", "Score", "Srand", "Enrichment", "pPvalue", "pQvalue", "OffsetFromStart")
  temp %<>% mutate(Length = END-START, pValue = 10^(-pPvalue))
  
  #filter out peaks with pValue > 10^-7
  #temp %<>% filter(pValue < 10^(-7))
  temp %>% arrange(CHR, START)
}, mc.cores = detectCores()) 


InputPeakAllCalled$GroupStat <- data.frame(Condition = names(InputPeakAllCalled$NarrowPeak),
                                                      TotalPeaksInGroup = lapply(InputPeakAllCalled$NarrowPeak, function(grp){
                                                        grp %>% .[!grepl("GL|hs", grp$CHR),] %>% nrow
                                                      }) %>% unlist,
                                                      TotalCoverage = lapply(InputPeakAllCalled$NarrowPeak, function(grp){
                                                        grp %>% .[!grepl("GL|hs", grp$CHR),] %>% .$Length %>% sum
                                                      }) %>% unlist)


InputPeakAllCalled$GroupStat %<>% mutate(TotalCovPercent = signif(100*TotalCoverage/(2.7*10^9), digits = 3))
                                           

############################# HTseq counts ######################################################################
HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub(".data.parkome.chipseq.results.bamfiles.0?", "", x)
  strsplit(x,"_")[[1]][1]
})
HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", paste0("X", as.character(Metadata$activemotif_id))))

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

#Filter out peaks with pValue > 10^-7 and peaks mapped to mitochondrial DNA/cotig regions
HTseqCounts %<>% filter(Geneid %in% InputPeakAllCalled$NarrowPeak$All$PeakName) %>% .[!grepl("GL|hs|MT", .$CHR),] %>% droplevels()

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("^X")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  

AllCalledData <- GetCountMatrixHTseq(HTseqCounts, OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29")


countMatrixPromotersAllCalled  <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot, collapseBy = "GeneAnnoType",CorMethod = "pearson",
                                                     FilterBy = "promoter", meta = AllCalledData$SampleInfo,
                                                     title = paste0("Sample correlation (promoter), Peaks - called together, ", Cohort))

pdf(paste0(ResultsPath, "SampleCorPromoter", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixPromotersAllCalled$HeatMap
dev.off()

StatAge <- sapply(levels(countMatrixPromotersAllCalled$Metadata$condition), function(Group){
  data = countMatrixPromotersAllCalled$Metadata %>% filter(condition == Group, !activemotif_id %in% c("57","39")) %>% droplevels()
  lm(RiP_NormMeanRatioOrg~age, data = data)
}, simplify = F)

ggplot(countMatrixPromotersAllCalled$Metadata %>% filter(!activemotif_id %in% c("57","39")), aes(age, RiP_NormMeanRatioOrg, color = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  geom_smooth(method = "lm", aes(color = condition, fill = condition), alpha = 0.3, size = 0.2) +
  geom_point() +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~condition, scales = "free_x")

StatNeuronProp <- sapply(levels(countMatrixPromotersAllCalled$Metadata$condition), function(Group){
  data = countMatrixPromotersAllCalled$Metadata %>% filter(condition == Group, !activemotif_id %in% c("57","39")) %>% droplevels()
  lm(RiP_NormMeanRatioOrg~NeuronProp, data = data)
}, simplify = F)

ggplot(countMatrixPromotersAllCalled$Metadata %>% filter(!activemotif_id %in% c("57","39")), aes(NeuronProp, RiP_NormMeanRatioOrg, color = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  geom_smooth(method = "lm", aes(color = condition, fill = condition), alpha = 0.3, size = 0.2) +
  geom_point() +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~condition, scales = "free_x")

ggplot(countMatrixPromotersAllCalled$Metadata %>% filter(!activemotif_id %in% c("57","39")), aes(age, NeuronProp,  color = condition)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_smooth(method = "lm", aes(color = condition, fill = condition), alpha = 0.3, size = 0.2) +
  geom_point() +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~condition, scales = "free_x")

ggplot(countMatrixPromotersAllCalled$Metadata %>% filter(!activemotif_id %in% c("57","39")), aes(NeuronProp, H3K27gapdh_Norm, color = condition)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_smooth(method = "lm", aes(color = condition, fill = condition), alpha = 0.3, size = 0.2) +
  geom_point() +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~condition, scales = "free_x")

ggplot(countMatrixPromotersAllCalled$Metadata %>% filter(!activemotif_id %in% c("57","39")), aes(age, H3K27gapdh_Norm, color = condition)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_smooth(method = "lm", aes(color = condition, fill = condition), alpha = 0.3, size = 0.2) +
  geom_point() +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  facet_wrap(~condition, scales = "free_x")

#Detect outliers
MedianCor <- apply(countMatrixPromotersAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]

Model = as.formula(" ~ sex + age + batch + NeuronProp")

DESeqOutAll_promoters <- RunDESeq(data = countMatrixPromotersAllCalled$countMatrix, UseModelMatrix = T, sampleToFilter = paste(names(Outlier), collapse = "|"),
                                  meta = countMatrixPromotersAllCalled$Metadata, normFactor = "MeanRatioAll",
                                  FullModel = Model, test = "Wald", FitType = "local")

                                    
DESegResultsSex_promotersAll <- GetDESeqResults(DESeqOutAll_promoters, coef = "sexM") %>% mutate(symbol = sapply(.$PeakName, function(x) {strsplit(x, "_")[[1]][1]}))
DESegResultsAge_promotersAll <- GetDESeqResults(DESeqOutAll_promoters, coef = "age") %>% mutate(symbol = sapply(.$PeakName, function(x) {strsplit(x, "_")[[1]][1]}))


###########################################################################################################
############ RERUN USING PEAKS BASED ON ALL THE SAMPLES, without collapsing peaks #########################
###########################################################################################################

countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",CorMethod = "pearson",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))

pdf(paste0(ResultsPath, "SampleCorAllPeaks", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixFullAllCalled$HeatMap
dev.off()

countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|^X"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("^X")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)


Model = as.formula(" ~ sex + age + batch + NeuronProp")

countMatrixDF <- AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)) %>% data.frame %>% select(matches("Peak|^X"))
countMatrixDF$MedianCount <- apply(countMatrixDF %>% select(matches("^X")), 1, median)
countMatrixDF %<>% mutate(NormCount = 200*MedianCount/Peak.width)
countMatrixDF$baseMean <- apply(countMatrixDF %>% select(matches("^X")), 1, mean)

countMatrix_filtered <- countMatrixDF %>% filter(NormCount > 5) %>% select(matches("^X")) %>% as.matrix()
rownames(countMatrix_filtered) <- as.character(countMatrixDF %>% filter(NormCount > 5) %>% .$PeakName)

################# Check the best normalization method ###############################

## Calculate normalization factors using RLE
DEStemp <- DESeqDataSetFromMatrix(countData = countMatrixDF %>% data.frame %>% select(matches("^X")) %>% as.matrix, colData = countMatrixFullAllCalled$Metadata, design = Model)
DEStemp <-  estimateSizeFactors(DEStemp)
DEStemp$RiP_RLE <- DEStemp$TotalCount/DEStemp$sizeFactor

MeltedDataAll <- DEStemp@colData %>% data.frame %>% gather(key = "MeasureType", value = "Value", TotalCount, library_size, RiP_NormAllCount, RiP_NormBackground, RiP_RLE, RiP_NormMeanRatioOrg, RiP_NormMeanRatioAll)
MeltedDataAll$MeasureType <- factor(MeltedDataAll$MeasureType, levels = c("TotalCount", "library_size",  "RiP_NormAllCount", "RiP_NormBackground", "RiP_RLE", "RiP_NormMeanRatioOrg", "RiP_NormMeanRatioAll"))
levels(MeltedDataAll$MeasureType) <- c("RiP", "LibrarySize",  "RiP/LibrarySize", "RiP/RoP", "RiP/sizeFactor", "RiP/MeanRatio", "RiP/MeanRatio2")

#Add individual vallues
xLabFun <- function(x) signif(as.numeric(as.character(x)), digits = 2)

StatNormMethod <- sapply(levels(MeltedDataAll$MeasureType), function(Type){
  data = MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm), activemotif_id != 57, MeasureType == Type) %>% droplevels()
  lm(H3K27gapdh_Norm~Value, data = data)
}, simplify = F)

MeltedDataAll %<>% mutate(MeasureType2 = MeasureType)

levels(MeltedDataAll$MeasureType2) <- sapply(levels(MeltedDataAll$MeasureType2), function(Type){
  temp <- StatNormMethod[[Type]] %>% summary()
  paste0(Type, " (r=", signif(temp$r.squared^0.5, digits = 2), ", p=",  signif(temp$coefficients[2,4], digits = 2), ")")
})


ggplot(MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm), activemotif_id != 57), aes(H3K27gapdh_Norm, Value, color = condition)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Counts", x = "WB, H3K27gapdh_normalized (Final)", title = "All samples Parkome") +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") + 
  #geom_smooth(method = "lm", color = "black") +
  scale_y_continuous(labels = xLabFun) +
  facet_wrap(~MeasureType2, scales = "free_y")
ggsave(paste0("WB_ChipSeqCor_FinalnormalizedToReplicates", Cohort, ".png"))



ggplot(MeltedDataAll %>% filter(MeasureType == "RiP/MeanRatio", batch != "H"), aes(age, Value, fill = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "Age", y = "Normalized RiP", title = "All samples") +
  geom_smooth(method = "lm", aes(fill = condition, color = condition), alpha = 0.3, size = 0.2) +
  geom_point(size = 2, aes(color = condition)) +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_y_continuous(labels = xLabFun) +
  facet_wrap(~condition, scales = "free_x")
ggsave(paste0("AgeRipCorrelation", Cohort, ".png"))

ggplot(MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm), MeasureType == "RiP/MeanRatio"), aes(age, H3K27gapdh_Norm, fill = condition)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "Age", y = "WB H3K27ac", title = "All samples") +
  geom_smooth(method = "lm", aes(fill = condition, color = condition), alpha = 0.3, size = 0.2) +
  geom_point(size = 2, aes(color = condition)) +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_y_continuous(labels = xLabFun) +
  facet_wrap(~condition, scales = "free_x")
ggsave(paste0("AgeWBCorrelation", Cohort, ".png"))


StatAgeGroup <- sapply(levels(MeltedDataAll$condition), function(group){
  data = MeltedDataAll %>% filter(!is.na(H3K27gapdh_Norm), condition == group) %>% droplevels()
  lm(age~H3K27gapdh_Norm, data = data)
}, simplify = F)


#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]

DESeqOutAll_Full <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = paste(names(Outlier), collapse = "|"),
                             FullModel = Model, test = "Wald", FitType = "local")
                                  

DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

#Repeat without PMI
Model2 = as.formula(" ~ condition + sex + age + batch")
DESeqOutAll_Full_noPMI <- RunDESeq(data = countMatrix_filtered, UseModelMatrix = T, 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = paste(names(Outlier), collapse = "|"),
                             FullModel = Model2, test = "Wald", FitType = "local")


DESegResultsSex_FullAll_noPMI <- GetDESeqResults(DESeqOutAll_Full_noPMI, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_noPMI <- GetDESeqResults(DESeqOutAll_Full_noPMI, coef = "age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_noPMI <- GetDESeqResults(DESeqOutAll_Full_noPMI, coef = "conditionPD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
