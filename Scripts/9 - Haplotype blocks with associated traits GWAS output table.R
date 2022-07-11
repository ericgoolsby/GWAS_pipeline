##### wrangling of the GWAS summary output files
library(tidyverse)
#source("functions.r")



colocate<-read.table("Tables/Colocate/Blocks/colocate_table.txt")
genecount<-read.csv("Tables/Colocate/Genes/genecount.csv")[,c(2,4)]
genelist<-read.csv("Tables/Colocate/Genes/Global_genelist.csv")
  genelist$X<-NULL


############## function that helps with settting a region as significant if it's only significant for one of it's component genome haplotype blocks
sig.sug.fun<-function (x) {
  if (sum("significant"%in%x)>0) {y<-"significant"} 
  if (sum("significant"%in%x)==0) {y<-"suggestive"}
  return(y)
}
##################


##### condense to single entry per region (collapse genome blocks) 
colocate<-colocate %>% dplyr::group_by(region,trait_env) %>%
                      dplyr::summarize(trait=trait[1],  
                                       env=env[1],
                                       pvalue=factor(sig.sug.fun(pvalue)),
                                       chromosome=chromosome[1],
                                       beta.sign=factor(sign(mean(beta.sign))))

colocate<-merge(colocate,genecount,by.x="region",by.y="colocate.block",all=T)



colocate.count<- colocate %>% filter(trait_env!="Tolerance_logdiff",trait_env!="Tolerance_water") %>%
                              group_by(region) %>%
                              dplyr::summarize(sig.nr=length(unique(trait_env[pvalue=="significant"])),
                                        sug.nr=length(unique(trait_env[pvalue=="suggestive"])),
                                        gene.nr=n[1],
                                        sig.traits=paste(as.character(unique(trait_env[pvalue=="significant"])),collapse=" / "),
                                        sug.traits=paste(as.character(unique(trait_env[pvalue=="suggestive"])),collapse=" / "))

#writing is csv so the following is unecessary
#colocate.count$region<-sub("-","_",colocate.count$region) # otherwise excel will change region ID's to dates



genelist<-merge(genelist,colocate.count,by.x="colocate.block",by.y="region")


write.csv(genelist,"Tables/Colocate/Genes/genelist-extended.csv",row.names=F)


write.csv(colocate.count,"Tables/Colocate/Genes/colocateCount.csv")

#significant count 
trait.count <- colocate %>% group_by(trait,env) %>% filter(pvalue=="significant") %>% 
                            dplyr::summarize(count=length(region)) %>% 
                            spread(env,count)

write.csv(trait.count,"Tables/Colocate/Genes/trait-regioncount.csv")

#suggestive couont
trait.count.suggestive <- colocate %>% group_by(trait,env) %>% filter(pvalue=="suggestive") %>% 
  dplyr::summarize(sug.count=length(region)) %>% 
  spread(env,sug.count)











#create genelist per traits

#replace empty cells with NA

traitgenelist <- genelist %>% mutate_all(na_if,"")

#separate sig and sug trait lists
sugtrait.genelist <- traitgenelist[!is.na(traitgenelist$sug.traits),]
sigtrait.genelist <- traitgenelist[!is.na(traitgenelist$sig.traits),]

#separate by rows
sugtrait.genelist <- separate_rows(sugtrait.genelist,"sug.traits",sep = " / ",convert = T)
sigtrait.genelist <- separate_rows(sigtrait.genelist,"sig.traits",sep = " / ",convert = T)

#remove extra data, add significance lable. rename trait column then combine

sugtrait.genelist <- subset(sugtrait.genelist, select = -sig.traits)
sugtrait.genelist$pvalue <- "suggestive"
names(sugtrait.genelist)[17]<-"traits"



sigtrait.genelist <- subset(sigtrait.genelist, select = -sug.traits)
sigtrait.genelist$pvalue <- "significant"
names(sigtrait.genelist)[17]<-"traits"


Trait.sig.sug.genelist<- rbind(sigtrait.genelist,sugtrait.genelist)



#separate into list of data frames based on column 

Separated.Trait.sig.sug.genelist <- split( Trait.sig.sug.genelist , f = Trait.sig.sug.genelist$traits )


names(Separated.Trait.sig.sug.genelist)<-paste(names(Separated.Trait.sig.sug.genelist),"_colocategenelist",sep = "")


#write list of dataframes to csv 

sapply(names(Separated.Trait.sig.sug.genelist), 
       function (x) write.csv(Separated.Trait.sig.sug.genelist[[x]], file=paste("Tables/Colocate/Genes/ColocateGenes/",x, ".csv", sep=""), row.names = FALSE )   )

#save separated GO terms for later 

Separated.Trait.sig.sug.genelist[is.na(Separated.Trait.sig.sug.genelist)] <- ""

sapply(names(Separated.Trait.sig.sug.genelist), 
       function (x) write.table(Separated.Trait.sig.sug.genelist[[x]][,c("locus_tag","Ontology_term")], file=paste("Tables/Colocate/Genes/ColocateGO/GO_",x, ".txt", sep=""),col.names=FALSE, row.names = FALSE, sep="\t" )   )




