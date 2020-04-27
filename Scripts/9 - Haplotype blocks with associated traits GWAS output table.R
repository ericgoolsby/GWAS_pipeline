##### wrangling of the GWAS summary output files
library(tidyverse)
#source("functions.r")



colocate<-read.table("Tables/Blocks/colocate_table.txt")
genecount<-read.csv("Tables/Genes/genecount.csv")[,c(2,4)]
genelist<-read.csv("Tables/Genes/genelist.csv")
  genelist$X<-NULL


############## function that helps with settting a region as significant if it's only significant for one of it's component genome haplotype blocks
sig.sug.fun<-function (x) {
  if (sum("significant"%in%x)>0) {y<-"significant"} 
  if (sum("significant"%in%x)==0) {y<-"suggestive"}
  return(y)
}
##################

##### condense to single entry per region (collapse genome blocks) 
colocate<-colocate %>% group_by(region,trait_env) %>%
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

colocate.count$region<-sub("-","_",colocate.count$region) # otherwise excel will change region ID's to dates

genelist<-merge(genelist,colocate.count,by.x="colocate.block",by.y="region")

write.csv(genelist,"Tables/Genes/genelist-extended.csv",row.names=F)


write.excel(colocate.count)

trait.count <- colocate %>% group_by(trait,env) %>% filter(pvalue=="significant") %>% 
                            dplyr::summarize(count=length(region)) %>% 
                            spread(env,count)

write.csv(trait.count,"Tables/Genes/trait-regioncount.csv")

trait.count.suggestive <- colocate %>% group_by(trait,env) %>% filter(pvalue=="suggestive") %>% 
  dplyr::summarize(sug.count=length(region)) %>% 
  spread(env,sug.count)

