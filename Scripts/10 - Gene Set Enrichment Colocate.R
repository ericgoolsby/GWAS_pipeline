
 library("topGO")
 
 
 
 #create file list of all files in Colocate GO & Colocate Genes for global 
 GO_term_files<-list.files(path = "Tables/Colocate/Genes/ColocateGO/", pattern = "global.txt")
 
 for (j in 1:length(GO_term_files)) {
    

 #import GO data
 geneID2GO<-readMappings(file=paste("Tables/Colocate/Genes/ColocateGO/",GO_term_files[j],sep = ""))

 
 geneNames<-names(geneID2GO)

 #create factor of geneNames to include all genes as significant
 geneList<-factor(as.integer(1:length(geneNames)%in% 1:length(geneNames)))
 
 #add a fake 2nd level
 geneList<-factor(geneList,levels = c(levels(geneList), 2))
 names(geneList) <- geneNames
 
 GOTYPES<-c("BP","MF")
 
 
 for (i in 1:length(GOTYPES)) {
    #create GOobject
    GOdata <-
       new(
          "topGOdata",
          ontology =  GOTYPES[i],
          allGenes = geneList,
          annot = annFUN.gene2GO,
          gene2GO = geneID2GO
       )
    
    resultFisher <-
       runTest(GOdata, algorithm = "classic", statistic = "fisher")
    
    resultKS <-
       runTest(GOdata, algorithm = "classic", statistic = "ks")
    resultKS.elim <-
       runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    try(allRes <- GenTable(
       GOdata,
       classicFisher = resultFisher,
       classicKS = resultKS,
       elimKS = resultKS.elim,
       orderBy = "elimKS",
       ranksOf = "classicFisher",
       topNodes = 20
    ))
    
   write.csv(allRes, paste("Tables/Colocate/Genes/ColocateGOResults/Global_Colocate_",GOTYPES[i],
    ".csv",sep = ""))
    
    
    pdf(paste(
       "Plots/Colocalization/GO/Global_Colocate_GO",
       GOTYPES[i],
       ".pdf",
       sep = ""
    ))
    showSigOfNodes(GOdata,
                   score(resultKS.elim),
                   firstSigNodes = 10,
                   useInfo = 'all')
    dev.off()
 }
 
 }
 
 
 
 
 #create file list of all files in Colocate GO & Colocate Genes for global 
 GO_term_files<-list.files(path = "Tables/Colocate/Genes/ColocateGO/", pattern = ".txt")
 GO_term_files<-grep(GO_term_files, pattern='global.txt', invert=T, value=T)
 
 
 
 GO_factor_files<-list.files(path = "Tables/Colocate/Genes/ColocateGenes/", pattern = ".csv")
 GO_factor_files<-grep(GO_factor_files, pattern='global.csv', invert = T, value=T)
 
 
 #filenames
 GO_resultnames<- gsub("^([^_]*_[^_]*)_.*$", "\\1", GO_factor_files)
 
 for (j in 1:length(GO_term_files)) {
    
    try({
    #import GO data
    geneID2GO<-readMappings(file=paste("Tables/Colocate/Genes/ColocateGO/",GO_term_files[j],sep = ""))
    
    
    geneNames<-names(geneID2GO)
    
 #import significance level
 GOfactor<- read.csv(paste("Tables/Colocate/Genes/ColocateGenes/",GO_factor_files[j],sep=""))
 
 
 GOfactor<-GOfactor[,c("locus_tag","pvalue")]
 
 mySignificantGenes<-GOfactor$locus_tag[GOfactor$pvalue=="significant"]
 #create gene list factor
  geneList<-factor(as.integer(1:length(geneNames) %in% 1:length(mySignificantGenes)))
  names(geneList) <- geneNames
if(length(geneNames)==length(mySignificantGenes)){
   #add a fake 2nd level
   geneList<-factor(geneList,levels = c(levels(geneList), 2))
   names(geneList) <- geneNames
   
}
  
  
    
  
  
  GOTYPES<-c("BP","MF")
  
  
  for (i in 1:length(GOTYPES)) {
     #create GOobject
     GOdata <-
        new(
           "topGOdata",
           ontology =  GOTYPES[i],
           allGenes = geneList,
           annot = annFUN.gene2GO,
           gene2GO = geneID2GO
        )
     
     resultFisher <-
        runTest(GOdata, algorithm = "classic", statistic = "fisher")
     
     resultKS <-
        runTest(GOdata, algorithm = "classic", statistic = "ks")
     resultKS.elim <-
        runTest(GOdata, algorithm = "elim", statistic = "ks")
     
     allRes <- GenTable(
        GOdata,
        classicFisher = resultFisher,
        classicKS = resultKS,
        elimKS = resultKS.elim,
        orderBy = "elimKS",
        ranksOf = "classicFisher",
        topNodes = 20
     )
     
     write.csv(allRes, paste("Tables/Colocate/Genes/ColocateGOResults/GO_Colocate_",GO_resultnames[j],"_",GOTYPES[i],
                             ".csv",sep = ""))
     
     
     pdf(paste(
        "Plots/Colocalization/GO/GO_Colocate",GO_resultnames[j],"_",
        GOTYPES[i],
        ".pdf",
        sep = ""
     ))
     showSigOfNodes(GOdata,
                    score(resultKS.elim),
                    firstSigNodes = 10,
                    useInfo = 'all')
     dev.off()
  }
    })
 } 
   
  
   
   
   
   
   