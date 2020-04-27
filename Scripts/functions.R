process_data <- function(trait_filename = "katie_data.csv",env_dat_to_merge)
{
  dat <- read.csv(paste("data/",trait_filename,sep=""),stringsAsFactors = FALSE)
  colnames(dat)[1] <- "SAM"
  colnames(dat) <- gsub(pattern = "_",replacement = "",x = colnames(dat))
  traits <- colnames(dat)[-1]
 # envs <- c("Dry","logdiff","Wet")
 if(!missing(env_dat_to_merge)){
   envs <- c("Env1","logdiff","Env2")
 }else{
   envs <- "common"
 }
  
  
  
   
  if(is.numeric(dat$SAM))
  {
    dat$SAM <- paste("SAM",sapply(3-nchar(dat$SAM),function(X) paste("",rep("0",X),sep="",collapse="")),dat$SAM,sep="")
  }
  if(!any(grepl("SAM001",dat$SAM))) stop("dat must be formatted as integers (SAM lines) or of the format e.g. SAM001, SAM002, etc.")
  load("data/lines.RDat")
  
  dat <- dat[match(x = lines,table = dat$SAM),]
  dat$SAM <- lines
  
  if(!missing(env_dat_to_merge))
  {
    env_dat_to_merge <- read.csv(paste("data/",env_dat_to_merge,sep=""),stringsAsFactors = FALSE)
    traits <- c(traits,unique(sapply(strsplit(colnames(env_dat_to_merge)[-1],split = "_"),function(X) X[[1]][1])))
  }
  
  writeLines(text = traits,"traits_to_run.txt")
  writeLines(text = envs,"environments_to_run.txt")
  
  new_names <- c("SAM",kronecker(X = envs,Y = colnames(dat)[-1],FUN = function(Y,X) paste(X,Y,sep="_")))
  #dat <- cbind(dat,dat[,-1]*0 + rnorm(n = length(unlist(dat[,-1]))),dat[,-1]*0 + rnorm(n = length(unlist(dat[,-1]))))
  #dat <- cbind(dat,dat[,-1],dat[,-1])
  colnames(dat) <- new_names
  
  if(!missing(env_dat_to_merge))
  {
   
    #this needs to be fixed to make the log difference 
    logdiff_betweenenviroments <-   
    dat <- cbind(dat,dat[,-1]*0 + rnorm(n = length(unlist(dat[,-1]))),dat[,-1]*0 + rnorm(n = length(unlist(dat[,-1]))))
    
     dat <- cbind(dat,env_dat_to_merge[,-1])
  }
  
  
  write.csv(x = dat,file = paste("data/",gsub(pattern = ".csv",replacement = "",x = trait_filename),"_for_pipeline.csv",sep=""),row.names = FALSE)
  
  prefs <- readLines("Scripts/original### Preferences ###")
  phen <-grep("replace1",prefs)
  prefs[phen] <- gsub(pattern = "replace1",replacement = paste(gsub(pattern = ".csv",replacement = "",trait_filename),"_for_pipeline.csv",sep=""),x = prefs[phen])
  writeLines(text = prefs,con = "Scripts/### Preferences ###")
}











# process_data <- function(trait_filename = "katie_data.csv",env_dat_to_merge)
# {
#   dat <- read.csv(paste("data/",trait_filename,sep=""),stringsAsFactors = FALSE)
#   colnames(dat)[1] <- "SAM"
#   colnames(dat) <- gsub(pattern = "_",replacement = "",x = colnames(dat))
#   traits <- colnames(dat)[-1]
#   envs <- c("Dry","logdiff","Wet")
#   
#   if(is.numeric(dat$SAM))
#   {
#     dat$SAM <- paste("SAM",sapply(3-nchar(dat$SAM),function(X) paste("",rep("0",X),sep="",collapse="")),dat$SAM,sep="")
#   }
#   if(!any(grepl("SAM001",dat$SAM))) stop("dat must be formatted as integers (SAM lines) or of the format e.g. SAM001, SAM002, etc.")
#   load("data/lines.RDat")
#   
#   dat <- dat[match(x = lines,table = dat$SAM),]
#   dat$SAM <- lines
#   
#   if(!missing(env_dat_to_merge))
#   {
#     env_dat_to_merge <- read.csv(paste("data/",env_dat_to_merge,sep=""),stringsAsFactors = FALSE)
#     traits <- c(traits,unique(sapply(strsplit(colnames(env_dat_to_merge)[-1],split = "_"),function(X) X[[1]][1])))
#   }
#   
#   writeLines(text = traits,"traits_to_run.txt")
#   writeLines(text = envs,"environments_to_run.txt")
#   
#   new_names <- c("SAM",kronecker(X = envs,Y = colnames(dat)[-1],FUN = function(Y,X) paste(X,Y,sep="_")))
#   dat <- cbind(dat,dat[,-1]*0 + rnorm(n = length(unlist(dat[,-1]))),dat[,-1]*0 + rnorm(n = length(unlist(dat[,-1]))))
#   #dat <- cbind(dat,dat[,-1],dat[,-1])
#   colnames(dat) <- new_names
#   
#   if(!missing(env_dat_to_merge))
#   {
#     dat <- cbind(dat,env_dat_to_merge[,-1])
#   }
#   
#   
#   write.csv(x = dat,file = paste("data/",gsub(pattern = ".csv",replacement = "",x = trait_filename),"_for_pipeline.csv",sep=""),row.names = FALSE)
#   
#   prefs <- readLines("Scripts/original### Preferences ###")
#   phen <-grep("replace1",prefs)
#   prefs[phen] <- gsub(pattern = "replace1",replacement = paste(gsub(pattern = ".csv",replacement = "",trait_filename),"_for_pipeline.csv",sep=""),x = prefs[phen])
#   writeLines(text = prefs,con = "Scripts/### Preferences ###")
# }

set_threshold <- function(method)
{
  if(method==1)
  {
    ntests <- 20562
  } else if(method==2)
  {
    ntests <- 50
  } else
  {
    if(is.na(as.numeric(method)))
    {
      stop("******************************************************************************************************
       Please set the argument method to one of the following values:\n\t1) for standard Bonferroni; assumes ntests = 20562  (p-value threshold = .05/ntests)
        2) Sets threshold to 0.001 by artificially assuming 50 tests.
        
        Alternatively, you can enter a custom number of tests in place of method.
       ***********************************************************************************************************")
    }
  }
  prefs <- readLines("Scripts/original### Preferences ###")
  thresh <-grep("replace2",prefs)
  prefs[thresh] <- gsub(pattern = "replace2",replacement = ntests,x = prefs[thresh])
  original_prefs <- readLines("Scripts/### Preferences ###")
  prefs[-thresh] <- original_prefs[-thresh]
  writeLines(text = prefs,con = "Scripts/### Preferences ###")
  
}

