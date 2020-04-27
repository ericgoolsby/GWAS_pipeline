#args <- as.numeric(na.exclude(as.numeric(commandArgs())))
# edit these lines as needed
### BEGIN SECTION 1 ###
trait_filename <- "SAMHPLCRetentiontimes200PlotPeakAreas.csv"


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




getwd()
system("chmod -R 755 ~/Documents/DocumentsPro2/UCFDissertation/GWAS_pipeline") # grant permissions (to avoid access denied errors)
pvalue_cutoff <- 1 # only change this for debugging; 1 = Bonferroni = 1; 2 = "suggested" 0.001 threshold


source("Scripts/functions.R")


process_data(trait_filename = trait_filename) #,env_dat_to_merge = "drought_and_rishi.csv")
set_threshold(method = pvalue_cutoff)
##############################
####### END SECTION ##########
##############################


##############################
### PIPELINE SCRIPTS 2-8 #####
##############################

  cat("\nBeginning Script 2.\n")
  source("Scripts/2 - Make manhattan plots.R")
  cat("\nCompleted Script 2.\n")
  cat("\nBeginning Script 2b.\n")
  source("Scripts/2b (optional) - Make single environment manhattan plots.R")
  cat("\nCompleted Script 2b.\n")
  

  cat("\nBeginning Script 3.\n")
  source("Scripts/3 - Blocks and heatmaps.R")
  cat("\nCompleted Script 3.\n")

  
  
  cat("\nBeginning Script 4.\n")
  source("Scripts/4 - Trait GWAS colocalization figure.R")
  cat("\nCompleted Script 4.\n")
  cat("\nBeginning Script 4c.\n")
  source("Scripts/4c - (optional) Per chromosome - Trait GWAS colocalization figure.R")
  cat("\nCompleted Script 4c.\n")
  

  cat("\nBeginning Script 5.\n")
  source("Scripts/5 - List genes in regions.R")
  cat("\nCompleted Script 5.\n")
  

  cat("\nBeginning Script 6.\n")
  attempt <- try(source("Scripts/6 - Manhattan with blocks.R"))
  if(inherits(x = attempt,what = "try-error"))
  {
    cat("Script 6 failed. X11 is a common error message if running on the cluster (note: this script currently requires an interactive R session.")
  } else
  {
    cat("\nCompleted Script 6.\n")
  }
  

  cat("\nBeginning Script 7.\n")
  source("Scripts/7 - Export PVE per trait to table.R")
  cat("\nCompleted Script 7.\n")
  

 cat("\nBeginning Script 8.\n")
  source("Scripts/8 - Show blocks on haplotype map.R")
  cat("\nCompleted Script 8.\n")



### BEGIN APPENDIX 1 ###
if(FALSE) # change to TRUE to install all required packages
{
  install.packages(c("data.table", "qqman", "tidyverse", "RColorBrewer", "ggpubr", 
                     "grid", "ggrepel", "gridExtra", "cowplot", "wesanderson", "corrr", 
                     "Hmisc", "ggdendro", "urltools", "scales"))
}
### END APPENDIX 1 ###
