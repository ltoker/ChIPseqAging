library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

source("ProjectScripts/ProjectFunctions.R")

ResultsDiscovery <- readRDS(paste0(DiscoveryResults, "/DESegResultsAge.L_FullAll.Rds"))
Deseq2OutDiscovery <- readRDS(paste0(DiscoveryResults, "/DESeqOutAll_Full.Rds"))

MarziAnalysis <- readRDS("Results_Aging_v35lift37/OutputMarzi.Rds")


DiseaseLDblocks <- list(AD = readRDS("data/LDblocks/LDsnpADRanges.Rds"),
                        PD = readRDS("data/LDblocks/LDsnpPDRanges.Rds"),
                        SCZ = readRDS("data/LDblocks/LDsnpSCZRanges.Rds"),
                        ASD = readRDS("data/LDblocks/LDsnpASDRanges.Rds"),
                        MS = readRDS("data/LDblocks/LDsnpMSRanges.Rds"))



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
      PeakInfo %>% filter(#baseMean >= BaseMeanQuantiles[Strat-1], baseMean < WidthQuantiles[Strat],
                          width >= WidthQuantiles[Strat-1], width < WidthQuantiles[Strat]) %>% .$PeakName
    } else if(Strat == 1){
      PeakInfo %>% filter(#baseMean < WidthQuantiles[1],
                          width < WidthQuantiles[1])  %>% .$PeakName
    } else if(Strat == StratNum){
      PeakInfo %>% filter(#baseMean >= BaseMeanQuantiles[Strat-1],
                          width >= WidthQuantiles[Strat-1])  %>% .$PeakName
    }
  })


  names(StratGroups) <- paste0("Strat_", 1:StratNum)


  RandomPeaks <- as.list(1:10000)
  names(RandomPeaks) <- paste0("Rand_", 1:10000)

  RandomPeaks <- lapply(RandomPeaks, function(x){
    lapply(StratGroups, function(Strat){
      sample(Strat, nrow(Diseasepeaks)/StratNum, replace = F)
    }) %>% unlist
  })
  if(length(RandomPeaks) < nrow(Diseasepeaks)){
    Extra = StratGroups[[length(StratGroups)]]
    Extra <- Extra[!Extra %in% RandomPeaks]
    RandomPeaks <- c(RandomPeaks, sample(Extra,
                                         nrow(Diseasepeaks) - length(RandomPeaks),
                                         replace = F)) 
  }

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

set.seed(666)

DiseaseLDEnrich <- lapply(DiseaseLDblocks, function(disease){
  GetLDblockEnrich(disease)
})

temp <- data.frame(RandomP  = DiseaseLDEnrich$AD$RandomSignif %>% unlist())
ggplot(temp, aes(-log10(RandomP))) +
  theme_minimal() +
  geom_density() +
  geom_vline(xintercept = -log10(DiseaseLDEnrich$AD$Hypergeometric), color = "red")

LDdiseaseDF <- sapply(names(DiseaseLDEnrich), function(Dis){
  data <- DiseaseLDEnrich[[Dis]]
  data.frame(Disease = Dis, LDpeaksAll = nrow(data$Diseasepeaks),
             LDpeaksSignif = nrow(data$SignifDiseasePeaks),
             pHypergeometric = scientific(data$Hypergeometric, digits = 2))
}, simplify = F) %>% rbindlist()

write.table(LDdiseaseDF, paste0(ResultsPath, "LDenrichment.tsv"),
            sep = "\t", row.names = F, col.names = T)

AllDiseasePeaks <- sapply(names(DiseaseLDEnrich), function(disease){
  DiseaseLDEnrich[[disease]]$Diseasepeaks %>% mutate(Disease = disease)
}, simplify = FALSE) %>% rbindlist %>% data.frame()

AllDiscovery <- ResultsDiscovery %>% select(PeakName, log2FoldChange, pvalue, padj) %>%
  filter(!duplicated(PeakName))

AllDiscovery <- merge(AllDiscovery, AllDiseasePeaks %>% select(PeakName, Disease), by = "PeakName", all.x = T)
AllDiscovery$Disease[is.na(AllDiscovery$Disease)] <- "Neither"

AllDiscovery$Disease <- factor(AllDiscovery$Disease, levels = c("Neither",  "ASD", "SCZ", "MS", "PD", "AD"))

MedianNeither <- AllDiscovery %>% filter(Disease == "Neither") %>% .$pvalue %>% log10(.) %>% "*"(-1) %>% median

AllDiscovery %<>% mutate(Direction = "Hypoacetylate")
AllDiscovery$Direction[AllDiscovery$log2FoldChange > 0] <- "Hyperacetylated"

ggplot(AllDiscovery , aes(Disease, -log10(pvalue))) +
  theme_minimal() +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_hline(yintercept = MedianNeither, color = "red", linetype = "dashed")


ggplot(AllDiscovery, aes(-log10(pvalue))) +
  theme_minimal() +
  geom_density(aes(fill = Disease, color = Disease), alpha = 0.3)
ggsave(paste0(ResultsPath, "LDlogPdist.pdf"), device = "pdf", width = 6, height = 5)

AllDiscovery$Microglia <- "Not Microglia-enriched"
AllDiscovery$Microglia[AllDiscovery$PeakName %in% (ChIPResults %>%
                                                     filter(symbol %in% AllMicroglia$feature) %>% .$PeakName)] <- "Microglia-enriched"
ggplot(AllDiscovery, aes(-log10(pvalue))) +
  theme_minimal() +
  geom_density(aes(fill = Disease, color = Disease), alpha = 0.3) +
  geom_vline(xintercept = (AllDiscovery %>% arrange(pvalue) %>%
                             filter(padj < 0.05) %>%
                             .$pvalue %>% max %>%
                             log10(.)*(-1)), lty = "dashed") +
  facet_wrap(Direction~Microglia)


#Look at the hypergeometric test after exclusion of Microglia peaks
PeakInfoNM <- ChIPResults %>% select(PeakName, Peak.CHR, Peak.START, Peak.END,
                                     Peak.width, baseMean, log2FoldChange, pvalue, padj) %>%
  filter(!is.na(padj), PeakName %in% (AllDiscovery %>%
                                        filter(Microglia == "Not Microglia-enriched") %>% .$PeakName)) %>%
  filter(!duplicated(PeakName))

names(PeakInfoNM) <- sapply(names(PeakInfoNM), function(x) gsub("Peak\\.", "", x)) 


DiseasepeaksNM <- DiseaseLDEnrich$AD$Diseasepeaks %>% filter(PeakName %in% PeakInfoNM$PeakName)

TotalTestableDiseasePeaksNM <- PeakInfoNM %>% filter(PeakName %in% DiseasepeaksNM$PeakName)
SignifDiseasePeaksNM <- DiseasepeaksNM %>% filter(padj < 0.05) 

TotalSignifPeaksNM <- PeakInfoNM %>% filter(padj < 0.05)

HyperNoMicroglia <- phyper(nrow(SignifDiseasePeaksNM)-1, nrow(TotalTestableDiseasePeaksNM),
                           nrow(PeakInfoNM)-nrow(TotalTestableDiseasePeaksNM),
                           nrow(TotalSignifPeaksNM), lower.tail = F)
