### script to make table of PVE of traits
library(tidyverse)
library(data.table)

colocate<-read.table("Tables/Colocate/Blocks/colocate_table.txt")

region.PVE<-colocate %>% filter(pvalue=="significant") %>%
                            group_by(trait_env,region) %>% dplyr::summarise(PVE=max(PVE))

trait.PVE<-region.PVE %>% group_by(trait_env) %>% dplyr::summarise(PVE=sum(PVE), count=length(region)) %>% separate(trait_env,into=c("trait","env"),sep="_")

PVE<-trait.PVE %>% gather("metric","value", -c(trait,env)) %>% unite("spreader",env,metric,sep="_") %>% spread(spreader,value)

write.table(PVE,"Tables/Colocate/Blocks/PVE.txt")
