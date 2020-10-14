source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")

ResultsDiscovery <- readRDS("AgingResults/DESegResultsAge.L_FullAll.Rds")
Deseq2OutDiscovery <- readRDS("AgingResults/DESeqOutAll_Full.Rds")

MarziAnalysis <- readRDS("AgingResults/OutputMarzi.Rds")




DiseaseLDblocks <- list(PD = readRDS("data/LDsnpPDRanges.Rds"),
                        AD = readRDS("data/LDsnpADRanges.Rds"),
                        SCZ = readRDS("data/LDsnpSCZRanges.Rds"),
                        ASD = readRDS("data/LDsnpASDRanges.Rds"),
                        MS = readRDS("data/LDsnpMSRanges.Rds"))



GetLDblockEnrich <- function(LDblocks, DESeqResult = ResultsDiscovery, StratNum = 5){
  PeakInfo <- DESeqResult %>% select(PeakName, Peak.CHR, Peak.START, Peak.END, Peak.width, baseMean, log2FoldChange, pvalue, padj) %>% filter(!duplicated(PeakName), !is.na(padj))
  names(PeakInfo) <- sapply(names(PeakInfo), function(x) gsub("Peak\\.", "", x)) 
  
  LDoverlap <- findOverlaps(LDblocks, PeakInfo %>% as(., "GRanges"))
  Diseasepeaks <- PeakInfo[subjectHits(LDoverlap),]
  Diseasepeaks %<>% filter(!duplicated(PeakName))
  
  TotalTestableDiseasePeaks <- PeakInfo %>% filter(PeakName %in% Diseasepeaks$PeakName)
  SignifDiseasePeaks <- Diseasepeaks %>% filter(padj < 0.05) 
  
  TotalSignifPeaks <- PeakInfo %>% filter(padj < 0.05)
  
  Signif <- phyper(nrow(SignifDiseasePeaks)-1, nrow(TotalTestableDiseasePeaks),
                   nrow(PeakInfo)-nrow(TotalTestableDiseasePeaks),
                   nrow(TotalSignifPeaks), lower.tail = F)
  
  
  BaseMeanQuantiles <- quantile(Diseasepeaks$baseMean, seq(1/StratNum, 1-1/StratNum, length.out = StratNum-1))
  WidthQuantiles <- quantile(Diseasepeaks$width, seq(1/StratNum, 1-1/StratNum, length.out = StratNum-1))

  StratGroups <- sapply(1:StratNum, function(Strat){
    if(Strat > 1 & Strat < StratNum){
      PeakInfo %>% filter(baseMean >= BaseMeanQuantiles[Strat-1], baseMean < WidthQuantiles[Strat],
                          width >= WidthQuantiles[Strat-1], width < WidthQuantiles[Strat]) %>% .$PeakName
    } else if(Strat == 1){
      PeakInfo %>% filter(baseMean < WidthQuantiles[1],
                          width < WidthQuantiles[1])  %>% .$PeakName
    } else if(Strat == StratNum){
      PeakInfo %>% filter(baseMean >= BaseMeanQuantiles[Strat-1],
                          width >= WidthQuantiles[Strat-1])  %>% .$PeakName
    }
  })


  names(StratGroups) <- paste0("Strat_", 1:StratNum)




  RandomPeaks <- as.list(1:1000)
  names(RandomPeaks) <- paste0("Rand_", 1:1000)

  RandomPeaks <- lapply(RandomPeaks, function(x){
    lapply(StratGroups, function(Strat){
      sample(Strat, nrow(Diseasepeaks)/StratNum, replace = F)
    }) %>% unlist
  })

  RandomSignif <- lapply(RandomPeaks, function(Random){
    RandomSignif <- PeakInfo %>% filter(PeakName %in% Random, padj < 0.05)
    phyper(nrow(RandomSignif)-1, nrow(TotalTestableDiseasePeaks),
           nrow(PeakInfo)-nrow(TotalTestableDiseasePeaks),
           nrow(TotalSignifPeaks), lower.tail = F)
  })
  
  
  
  return(list(Diseasepeaks = Diseasepeaks,
              SignifDiseasePeaks = SignifDiseasePeaks,
              Hypergeometric = Signif,
              RandomSignif = RandomSignif,
              RandomPeaks = RandomPeaks))
}

seed(666)

DiseaseLDEnrich <- lapply(DiseaseLDblocks, function(disease){
  GetLDblockEnrich(disease)
})



AllDiseasePeaks <- sapply(names(DiseaseLDEnrich), function(disease){
  DiseaseLDEnrich[[disease]]$Diseasepeaks %>% mutate(Disease = disease)
}, simplify = FALSE) %>% rbindlist %>% data.frame()

AllDiscovery <- ResultsDiscovery %>% select(PeakName, log2FoldChange, pvalue, padj) %>% filter(!duplicated(PeakName), !is.na(padj))

AllDiscovery <- merge(AllDiscovery, AllDiseasePeaks %>% select(PeakName, Disease), by = "PeakName", all.x = T)
AllDiscovery$Disease[is.na(AllDiscovery$Disease)] <- "Neither"

AllDiscovery$Disease <- factor(AllDiscovery$Disease, levels = c("Neither", "MS", "ASD", "SCZ", "PD", "AD"))

MedianNeither <- AllDiscovery %>% filter(Disease == "Neither") %>% .$pvalue %>% log10(.) %>% "*"(-1) %>% median


ggplot(AllDiscovery , aes(Disease, -log10(pvalue))) +
  #geom_violin(width = 0.4) +
  geom_boxplot(width = 0.1) +
  geom_hline(yintercept = MedianNeither, color = "red", linetype = "dashed")


ADsignifPeakCounts <- sapply(SignifADsnpPeaks$PeakName, function(peak){
  temp <- plotCounts(Deseq2OutDiscovery, gene = peak, intgroup = "Age", normalized = T, transform = T, returnData = T)
  temp$Peak <- peak
  temp$Disease <- "AD"
  temp$count <- scale(temp$count)
  temp
}, simplify = F) %>% rbindlist() %>% data.frame()

ggplot(ADsignifPeakCounts, aes(Age, count, color = Peak))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_jitter() +
  geom_smooth(se = FALSE) + facet_wrap(~Peak)


PDsignifPeakCounts <- sapply(SignifPDsnpPeaks$PeakName, function(peak){
  temp <- plotCounts(Deseq2OutDiscovery, gene = peak, intgroup = "Age", normalized = T, transform = T, returnData = T)
  temp$Peak <- peak
  temp$Disease <- "PD"
  temp$count <- scale(temp$count)
  temp
}, simplify = F) %>% rbindlist() %>% data.frame()

ggplot(PDsignifPeakCounts, aes(Age, count, color = Peak))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_jitter() +
  geom_smooth(se = FALSE) + facet_wrap(~Peak)



SCZsignifPeakCounts <- sapply(SignifSCZsnpPeaks$PeakName, function(peak){
  temp <- plotCounts(Deseq2OutDiscovery, gene = peak, intgroup = "Age", normalized = T, transform = T, returnData = T)
  temp$Peak <- peak
  temp$Disease <- "SCZ"
  temp$count <- scale(temp$count)
  temp
}, simplify = F) %>% rbindlist() %>% data.frame()

ggplot(SCZsignifPeakCounts, aes(Age, count, color = Peak))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_jitter() +
  geom_smooth(se = FALSE) + facet_wrap(~Peak)






ggplot(DiseaseGenes, aes(Disease, -log10(padj))) +
  geom_violin() + geom_boxplot(outlier.shape = NA, width = 0.1)
