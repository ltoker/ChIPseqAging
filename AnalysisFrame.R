BiocManager::install("devtools", force = T)
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

source("ProjectScripts/ProjectFunctions.R")


##### Main ChIPseq analysis ############
Cohort = "Aging"

source("ProjectScripts/Annotations.R") #This script defines the  assembly version 
                                       #of the annotations, adjust accordingly!

ResultsPath = paste0("Results_", Cohort, "_", strsplit(AssemblyFilename, "\\.")[[1]][2])
source("ProjectScripts/AnalysisAging.R")

rm(list = ls(all.names = TRUE))
closeAllConnections()

#Run replication 
Cohort = "MarziAD_OurPeaks"
source("ProjectScripts/Annotations.R")

DiscoveryResults = paste0("Results_Aging_", strsplit(AssemblyFilename, "\\.")[[1]][2])
ResultsPath = paste0("Results_", Cohort, "_", strsplit(AssemblyFilename, "\\.")[[1]][2])

source("ProjectScripts/ReplicationMarzi.R")

rm(list = ls(all.names = TRUE))
closeAllConnections()

##### LD block analysis #########

#Creating LD blocks based on the disease-associated SNPs

source("ProjectScripts/GetLDblocks.R") #If not using the existing ASDsnps.txt,
                                       #this would require manual editing of the file after creation


#Running the LD enrichment analysis
Cohort = "Aging"

source("ProjectScripts/Annotations.R") 

ResultsPath = paste0("Results_", Cohort, "_", strsplit(AssemblyFilename, "\\.")[[1]][2] ,"/")

DiscoveryResults = paste0("Results_Aging_", strsplit(AssemblyFilename, "\\.")[[1]][2],"/")
source("ProjectScripts/LDenrichmentAnalysis.R")

rm(list = ls(all.names = TRUE))
closeAllConnections()

#####  RNAseq analysis #########
Cohort = "Aging"
source("ProjectScripts/Annotations.R")

ResultsPath = paste0("Results_", Cohort, "_",
                     strsplit(AssemblyFilename, "\\.")[[1]][2] ,"/RNAseq/")

DiscoveryResults = paste0("Results_Aging_", strsplit(AssemblyFilename, "\\.")[[1]][2], "/")
SampleIDcol = "RNAseq_id_ParkOme2"

source("ProjectScripts/RNAseqAnalysis.R")

#Reanalysis Galatro et al.

Cohort = "Galatro"

source("ProjectScripts/Annotations.R")

ResultsPath = paste0("Results_", Cohort, "_",
                     strsplit(AssemblyFilename, "\\.")[[1]][2])
SampleIDcol = "Filename" 
