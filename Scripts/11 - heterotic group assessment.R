install.packages("pastecs")
install.packages("ggExtra")
library(ggExtra)
library(pastecs)
library(ggplot2)

#read in the trait data that was used for GWAS
alltraitdata<-read.csv(paste0(getwd(),"/data/",trait_filename))                       
#set colnames -1 as a list 
Trait_names<-colnames(alltraitdata[,-1])
Trait_names
#get all files that have log.txt from assoc tables

GEMMA_Logs<-list.files(path = paste0(getwd(),"/Tables/Assoc_files"), pattern = "log.txt")

#set a mastersheet

i<-1
#get the name of the log files
triallog<-paste0(getwd(),"/Tables/Assoc_files/",GEMMA_Logs[1])

#import log data except for model betas
logdata<-read.table(triallog,skip = 16, nrows=12, comment.char = "", sep = "#")
Mastershet<- as.data.frame(c(gsub("_.*","",GEMMA_Logs[1]),as.numeric(gsub(".*= ","",logdata[,3]))))


for (i in 1:length(GEMMA_Logs)) {
#get the name of the log files
triallog<-paste0(getwd(),"/Tables/Assoc_files/",GEMMA_Logs[i])

#import log data except for model betas
logdata<-read.table(triallog,skip = 16, nrows=12, comment.char = "", sep = "#")
Mastershet[,i+1]<- c(gsub("_.*","",GEMMA_Logs[i]),as.numeric(gsub(".*= ","",logdata[,3])))
#Mastershet<-rbind(Mastersheet,logentry)

}
View(Mastershet)
Mastershet<-t(Mastershet)[-1,]
rownames(Mastershet)<- Mastershet[,1]
colnames(Mastershet)<-c("Trait","GenotypeObs","PhenotypeObs","Covariates","Phenotypes","GenotypeSNPs","PhenotypeSNPs","REMLENULL","MLENULL","PVENull","SEpve","Vkin","Vpop")



#get descriptive statistices for All traits
ugly<-t(stat.desc(alltraitdata)[,-1])

All_Descriptive_Statistics<- merge(Mastershet, ugly, by=0, all=TRUE)

All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)]<-sapply(All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)],as.character)
All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)]<-sapply(All_Descriptive_Statistics[,3:ncol(All_Descriptive_Statistics)],as.numeric)



write.csv(All_Descriptive_Statistics, paste0(getwd(),"/Tables/Descriptive_Statistics",trait_filename))

View(All_Descriptive_Statistics)
#read in trait metadata

TRAITMETA<-read.csv(paste0(getwd(),"/data/AllTraitsforGWAS_2020_METADATA.csv"))

View(TRAITMETA)
#create subsetdata

data<-data.frame(var1=All_Descriptive_Statistics$PVENull,META1=TRAITMETA$META1,META2=TRAITMETA$META2, META3=TRAITMETA$META3)
#create labels for graphs

bottomlabel<- expression(italic(h)^2~""[SNP])
 


# classic plot :

plot_multi_histogram <- function(df, feature, label_column,bottomlabel) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black",bins=50) +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=bottomlabel, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))
  plt + theme_classic()
}

plot_multi_histogram(data,"var1","META3",bottomlabel = bottomlabel)
?ggplot

#read in significant effect sizes

Trait.EffectSIZE.Signif<- read.table("Tables/Colocate/Blocks/traits_to_genomeblocks_signif.txt",header = T)


#add metadata column 

Trait.EffectSIZE.Signif.META<-merge(Trait.EffectSIZE.Signif,TRAITMETA, by="trait")

View(Trait.EffectSIZE.Signif.META)

plot_multi_histogram(Trait.EffectSIZE.Signif.META,"PVE","META3",bottomlabel = "Effect Size")


#read in suggested effect sizes


Trait.EffectSIZE.sugest<- read.table("Tables/Colocate/Blocks/traits_to_genomeblocks_sugest.txt",header = T)


#add metadata column 

Trait.EffectSIZE.sugest.META<-merge(Trait.EffectSIZE.sugest,TRAITMETA, by="trait")

View(Trait.EffectSIZE.sugest.META)

plot_multi_histogram(Trait.EffectSIZE.sugest.META,"PVE","META3",bottomlabel = "Effect Size")





#plot additive effect size against heritibility

Trait.EffectSIZE.Additive<- read.table("Tables/Colocate/Blocks/PVE.txt",header = T)


Trait.EffectSIZE.Additive.meta<-merge(Trait.EffectSIZE.Additive,All_Descriptive_Statistics, by.x = "trait", by.y="Trait")

Trait.EffectSIZE.Additive.meta<-merge(Trait.EffectSIZE.Additive.meta,TRAITMETA,by="trait")

View(Trait.EffectSIZE.Additive.meta)

plot(Trait.EffectSIZE.Additive.meta$common_PVE,Trait.EffectSIZE.Additive.meta$PVENull)




ggplot(Trait.EffectSIZE.Additive.meta, aes(x=common_PVE, y=PVENull,col=META3)) + geom_point() +theme_classic()





