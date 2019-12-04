FullMeta <- read.table("meta/Neuromics_Brain_bank_v3_new_IF.csv", header = T, sep = "\t", quote = '"')
Metadata <- read.table("meta/AgingChipMetadata.csv", header = T, sep = "\t", na.strings = "NA")
Metadata$FinalBatch <- Metadata$sequencing_2 %>% as.character()
Metadata$FinalBatch[is.na(Metadata$FinalBatch)] <- as.character(Metadata$sequencing_1)[is.na(Metadata$FinalBatch)]
Metadata$FinalBatch %<>% factor()

Metadata$FinalFileName <- Metadata$filename_2 %>% as.character()
Metadata$FinalFileName[is.na(Metadata$FinalFileName)] <- as.character(Metadata$filename_1)[is.na(Metadata$FinalFileName)]

Metadata$PeakFileSample <- sapply(Metadata$FinalFileName, function(x){
  gsub("fastq", "narrowPeak", x)
}) %>% as.character()
Metadata$PeakFileAll <- "all.narrowPeak.gz"

Metadata$SampleID <- sapply(Metadata$FinalFileName, function(x){
  paste0("X", strsplit(as.character(x), "_")[[1]][1]) %>% gsub("X0", "X", .)
})


#Filtering out sample 56, since there was not enough material to resequence it an the Down syndrome infant
Metadata %<>% filter(!SampleID  %in% c("X56", "X72")) %>% droplevels()

#Add the biobank metadata
Metadata <- merge(Metadata, FullMeta %>% select(Biobank.ID, Cohort, Age, Sex, PMI, PH, Clinical.diagnosis, Pathology.diagnosis, braak.tau, amyloid, braak.LB, apoE, Cause.of.death, Other.diseases, Drugs.last.year.of.life), by.x = "BankID", by.y = "Biobank.ID", sort = F)
Metadata$Age <- as.numeric(as.character(Metadata$Age))
Metadata %<>% mutate(Agef = cut(Age, 5, ordered_result = T))
Metadata$AgeGroup <- sapply(Metadata$Age, function(x){
  if(x < 20){
    "Young"
  } else if(x > 50){
    "Old"
  } else {
    "Middle"
  }
}) %>% factor(levels = c("Middle", "Young", "Old"))

Metadata$Cohort2 <- sapply(as.character(Metadata$Cohort), function(x){
  if(x == "Netherlands Brain Bank"){
    "NBB"
  } else if(x == "Neuromics Tissue Bank"){
    "PV"
  } else {
    NA
  }
})


QC <- read.table("meta/all_samples_cc.txt", header = T, sep = "\t")
QC %<>% mutate(FinalFileName = Filename)
QC$FinalFileName <- sapply(QC$FinalFileName, function(x){
  gsub("tagAlign", "fastq", x)
})

Metadata <- merge(Metadata, QC %>% select(numReads, estFragLen, QualityTag, FinalFileName), by = "FinalFileName", sort = F)


#ADD just for now the additional sample 60, run originally in batch 4
# temp <-  Metadata %>% filter(sample_id == "60")
# temp %<>% mutate(FinalBatch = sequencing_1,
#                  FinalFileName = filename_1,
#                  PeakFileSample = gsub("fastq", "narrowPeak", FinalFileName),
#                  SampleID = paste0(SampleID, "batch4"))
# Metadata <- rbind(Metadata, temp)

