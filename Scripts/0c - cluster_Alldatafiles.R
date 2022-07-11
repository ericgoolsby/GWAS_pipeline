args <- as.numeric(na.exclude(as.numeric(commandArgs())))
# edit these lines as needed
### BEGIN SECTION 1 ###



#lapply(c("data.table", "qqman", "tidyverse", "RColorBrewer", "ggpubr", "grid", "ggrepel", "gridExtra", "cowplot", "wesanderson", "corrr", "dplyr", "Hmisc", "ggdendro", "urltools", "scales"),library,character.only=TRUE)

### END SECTION 1 ###

##############################
### RUN BUT DO NOT EDIT ######
##### THIS SECTION ###########
##############################
##### BEGIN SECTION ##########
##############################

requiredPackages <- c("data.table", "qqman", "tidyverse", "RColorBrewer", "ggpubr", "grid", "ggrepel", "gridExtra", "cowplot", "wesanderson", "corrr", "dplyr", "Hmisc", "ggdendro", "urltools", "scales")
for(Packagesneeded in requiredPackages){
   if(!require(Packagesneeded,character.only = TRUE)) install.packages(Packagesneeded)
   library(Packagesneeded,character.only = TRUE)
}


BIOPackages <- c("topGO")
for(Packagesneeded in BIOPackages){
   if(!require(Packagesneeded,character.only = TRUE)) BiocManager::install("topGO")
   library(Packagesneeded,character.only = TRUE)
}

system("chmod -R 755 ../GWASPIPELINE") # grant permissions (to avoid access denied errors)
pvalue_cutoff <- 1 # only change this for debugging; 1 = Bonferroni = 1; 2 = "suggested" 0.001 threshold
source("Scripts/functions.R")



datafiles<-list.files(path = "data/", pattern=".csv")
datafiles<-datafiles[!grepl("pipeline", datafiles)]

print(datafiles[1])
datafilenames<- sub(".(csv)", "", datafiles)
LL<-2
for (LL in 1:length(datafiles)) {
   trait_filename<-datafiles[LL]
   environmentalfile<-list.files(path = "data/Environment2", pattern=substr(trait_filename,start=1,stop = 5))
   #if a matching file does not exist proceed treating this as a common garden experiment 
   if((is.na(environmentalfile[1]))){
      process_data(trait_filename = trait_filename) 
   }else{
      process_data(trait_filename = trait_filename, env_dat_to_merge = environmentalfile[1]) 
  }
   
   
   
   
   
   #,env_dat_to_merge = "drought_and_rishi.csv")
   set_threshold(method = pvalue_cutoff)
   ##############################
   ####### END SECTItrait.data<-pheno.data[, select_cols]N ##########
   ##############################
   
   
   
   
   
   ##############################
   ### PIPELINE SCRIPTS 1-8 #####
   ##############################
   if(any(args==1))
   {
      cat("\nBeginning Script 1.\n")
      system("rm Tables/Assoc_files/*")
      source("Scripts/1 - Phenotype to GEMMA.R")
      cat("\nCompleted Script 1.\n")
   }
   if(any(args==2))
   {
      cat("\nBeginning Script 2.\n")
      source("Scripts/2 - Make manhattan plots.R")
      cat("\nCompleted Script 2.\n")
      cat("\nBeginning Script 2b single manhattan plots.\n")
      source("Scripts/2b (optional) - Make single environment manhattan plots.R")
      cat("\nCompleted Script 2b.\n")
      
   }
   if(any(args==3))
   {
      cat("\nBeginning Script 3.\n")
      source("Scripts/3 - Blocks and heatmaps.R")
      cat("\nCompleted Script 3.\n")
      cat("\nBeginning Script Single Trait 3.\n")
      source("Scripts/3b (SINGLETRAIT) SNPS in blocks.R")
      cat("\nCompleted Script 3.\n")
   }
   if(any(args==4))
   {
      cat("\nBeginning Script 4.\n")
      source("Scripts/4 - Trait GWAS colocalization figure.R")
      cat("\nCompleted Script 4.\n")
      cat("\nBeginning Script 4c.\n")
      source("Scripts/4c - (optional) Per chromosome - Trait GWAS colocalization figure.R")
      cat("\nCompleted Script 4c.\n")
      
   }
   if(any(args==5))
   {
      cat("\nBeginning Script 5.\n")
      source("Scripts/5 - List genes in regions.R")
      cat("\nCompleted Script 5.\n")
      
   }
   if(any(args==6))
   {
      cat("\nBeginning Script 6.\n")
      attempt <- try(source("Scripts/6 - Manhattan with blocks.R"))
      if(inherits(x = attempt,what = "try-error"))
      {
         cat("Script 6 failed. X11 is a common error message if running on the cluster (note: this script currently requires an interactive R session.")
      } else
      {
         cat("\nCompleted Script 6.\n")
      }
      
   }
   if(any(args==7))
   {
      cat("\nBeginning Script 7.\n")
      source("Scripts/7 - Export PVE per trait to table.R")
      cat("\nCompleted Script 7.\n")
      
   }
   if(any(args==8))
   {
      cat("\nBeginning Script 8.\n")
      source("Scripts/8 - Show blocks on haplotype map.R")
      cat("\nCompleted Script 8.\n")
   }
   if(any(args==9))
   {
      cat("\nBeginning Script 9.\n")
      source("Scripts/9 - Haplotype blocks with associated traits GWAS output table.R")
      cat("\nCompleted Script 9.\n")
   }
   
   if(any(args==10))
   {
      cat("\nBeginning Script 10.\n")
      source("Scripts/10 - Gene Set Enrichment Colocate.R")
      source("Scripts/10a - Gene Set Enrichment Colocate.R")
      cat("\nCompleted Script 10.\n")
   }
   
   ### BEGIN APPENDIX 1 ###
   if(FALSE) # change to TRUE to install all required packages
   {
      install.packages(c("data.table", "qqman", "tidyverse", "RColorBrewer", "ggpubr", 
                         "grid", "ggrepel", "gridExtra", "cowplot", "wesanderson", "corrr", 
                         "Hmisc", "ggdendro", "urltools", "scales"))
   }
   ### END APPENDIX 1 ###
   
   
   
   file.rename("Plots",paste(datafilenames[LL],"_Plots",sep=""))
   file.rename("Tables",paste(datafilenames[LL],"_Tables",sep=""))
   
   }
   