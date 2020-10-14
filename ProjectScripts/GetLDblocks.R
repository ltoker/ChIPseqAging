source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
#LDlinkR access token: 6c4c60024204
Token = "6c4c60024204"
BiocManager::install("LDlinkR")
library(LDlinkR)
packageF("GenomicRanges")

packageF("tabulizer")



#Add SNPs from releveant publications
GetR2 <- function(SNPs){
  SNPs$CHR <- paste0("chr", SNPs$CHR)
  
  snpR2 <- sapply(as.character(SNPs$SNP), function(x){
    LDproxy(x, "EUR", "r2",
            token = Token) 
  }, simplify = F)
  return(snpR2)
}

GetLDblock <- function(snpR2, SNPs, threshold = 0.6){
  LDsnp <- lapply(snpR2, function(x){
    if(nrow(x) > 1){
      temp <- x %>% filter(R2 >= threshold)
      temp$Chr <- sapply(temp$Coord, function(i) strsplit(as.character(i), ":")[[1]][1]) 
      temp$Loc <- sapply(temp$Coord, function(i) strsplit(as.character(i), ":")[[1]][2]) %>% as.numeric()
      temp
    } else {
      NA
    }
  })
  
  LDsnpRanges <- sapply(names(LDsnp), function(snp){
    if(length(LDsnp[[snp]]) > 1){
      Data = LDsnp[[snp]]
      Start = Data %>% arrange(Distance) %>% select(RS_Number, Chr, Loc, Distance, R2) %>% .[1,]
      End = Data %>% arrange(desc(Distance)) %>% select(RS_Number, Chr, Loc, Distance, R2) %>% .[1,]
      data.frame(SNPbait = snp, SNPstart = Start$RS_Number, SNPend = End$RS_Number,
                 SNPstartR2 = Start$R2, SNPendR2 = End$R2,
                 CHR = unique(Data$Chr), START = Start$Loc, END = End$Loc,
                 Width = End$Loc - Start$Loc)
    } else {
      Data = SNPs %>% filter(SNP == snp)
      data.frame(SNPbait = snp, SNPstart = NA, SNPend = NA,
                 SNPstartR2 = NA, SNPendR2 = NA,
                 CHR = Data$CHR, START = Data$BP, END = Data$BP,
                 Width = NA)
    }
  }, simplify = F) %>% rbindlist %>% data.frame() %>% droplevels() %>% as(., "GRanges")
}


ADsnps <- read.table("data/ADsnps.txt", header = T, sep = "\t")
ADsnpR2 <- GetR2(ADsnps)
LDsnpADRanges <- GetLDblock(ADsnpR2, ADsnps)


saveRDS(LDsnpADRanges, "data/LDsnpADRanges.Rds")
saveRDS(ADsnpR2, "data/LDsnpAD.Rds")



PDsnps <- read.table("data/PD_snps.txt", header = T, sep = "\t")
PDsnpR2 <- GetR2(PDsnps)
LDsnpPDRanges <- GetLDblock(PDsnpR2, PDsnps)



saveRDS(LDsnpPDRanges, "data/LDsnpPDRanges.Rds")
saveRDS(PDsnpR2, "data/LDsnpPD.Rds")


SCZsnps <- read.table("data/SCZsnps.txt", header = T, sep = "\t")

#Filter the "official" SNPs
SCZsnps <- SCZsnps[grepl("rs", SCZsnps$SNP),]

SCZsnpR2 <- GetR2(SCZsnps)

#Get the SNP coordinates since they were not provided in the original publication
SCZsnps$BP <- sapply(names(SCZsnpR2), function(SNP){
  if(SNP != "rs7907645" ){
    if(length(SCZsnpR2[[SNP]]) > 1){
      if(SNP %in% SCZsnpR2[[SNP]]$RS_Number){
        temp <- SCZsnpR2[[SNP]] %>% filter(RS_Number == SNP) %>% .$Coord %>% as.character()
      } else {
        temp <- SCZsnpR2[[SNP]] %>% filter(R2 == 1) %>% .$Coord %>% as.character()
      }
      strsplit(temp, ":")[[1]][2] %>% as.numeric()
    } else {
      NA
    }
  } else {
    104423800
  }
}) 

LDsnpSCZRanges <- GetLDblock(SCZsnpR2, SCZsnps)
LDsnpSCZRanges %>% head

saveRDS(LDsnpSCZRanges, "data/LDsnpSCZRanges.Rds")
saveRDS(SCZsnpR2, "data/LDsnpSCZ.Rds")


ASDsnps <- extract_tables("data/Groove_ADsup.pdf", pages = 45, method = "stream", output = "data.frame") %>%
  rbindlist() %>% data.frame()

ASDsnps <- rbind(ASDsnps, extract_tables("data/Groove_ADsup.pdf", pages = 46, method = "stream", output = "data.frame") %>%
  rbindlist() %>% data.frame()) %>% data.frame()

ASDsnps %<>% select(SNP, CHR, BP)
ASDsnps %<>% mutate(CHR = paste0("chr", CHR))
ASDsnps$BP <- sapply(ASDsnps$BP, function(x){
  gsub(",", "", x) %>% as.numeric()
})

write.table(ASDsnps, "data/ASDsnps.txt", sep = "\t", col.names = T, row.names = F)

ASDsnps <- read.table("data/ASDsnps.txt", sep = "\t", header = T)

ASDnpR2 <- GetR2(ASDsnps)

LDsnpASDRanges <- GetLDblock(ASDnpR2, ASDsnps)

saveRDS(LDsnpASDRanges, "data/LDsnpASDRanges.Rds")
saveRDS(ASDnpR2, "data/LDsnpASD.Rds")

#MSP - Supp table S7 Patsopoulos et al. PMID  31604244
MSsnps <- read.table("data/MSsnps.txt", sep = "\t", header = T) 
NSsnpR2 <- GetR2(MSsnps)

LDsnpMSRanges <- GetLDblock(MSsnpR2, MSsnps)

saveRDS(LDsnpMSRanges, "data/LDsnpMSRanges.Rds")
saveRDS(MSsnpR2, "data/LDsnpMS.Rds")
