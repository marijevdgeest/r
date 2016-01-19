##########################################################
#Authors: Jennifer Di Tomasso & Marije van der Geest     #
#                                                        #
#Performs cell counts/proportions prediction             #
#                                                        #
##########################################################

########
# Libs #
########
library(plyr)
library(foreach, lib.loc="~/r_libs/")
library(glmnet, lib.loc="~/r_libs/")


#Create connection
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]","[PASSWORD]")


#Function: cell counts pedictor 
pred.cc <- function(expression){

  expression <- as.matrix(expression)
  cell.type <- c("Neutrophils","Lymphocytes","Monocytes","Eosinophils","Basophils","Neutrophils%","Lymphocytes%","Monocytes%","Eosinophils%","Basophils%")
  prediction <- matrix(nrow=ncol(expression),ncol=10,dimnames=list(colnames(expression),cell.type))
  
  #For each cell type
  for (ct in 1:10){
    pred <- NULL
    #For each model per celltype
    for (m in 1:100){
      #Get model
      load(paste0("/srv/molgenis/.molgenis/omx/data/models/out_",cell.type[ct],"_",m,".Rdata")) ####check pathway
      pred.counts <- predict(out$fit1, t(expression), s=out$cvfit$lambda.min)
      pred <- cbind(pred,pred.counts)
    }
    prediction[,ct] <- apply(pred, 1, mean)
  }
  return(prediction)
  
}

sampleGenes <- molgenis.get("[SAMPLEGENES]", num=100000)

demo <- data.frame(sampleGenes, row.names=1)
upload <- molgenis.get("${importedEntity}", num=100000)
upload.fixed <- data.frame(upload, row.names=1)
#list all genes
gene.list <- rownames(demo) 
#filter out genes from example dataset in order to list the missing genes
gene.diff <- setdiff(gene.list,rownames(upload.fixed)) 

#Create matrix for output
mat.zero <- matrix(0,ncol=ncol(upload.fixed),nrow=(length(gene.list)-nrow(upload.fixed)),dimnames = list(gene.diff,colnames(upload.fixed))) #create zero matrix with all labels
#create final matrix for the predictor
final <- rbind(upload.fixed, mat.zero) 
#re-order gene names
final.ordered <- final[gene.list,] 


#Start prediction
output <- pred.cc(as.matrix(final.ordered))

#Format output
output2 <- cbind(rownames(output), output)
rownames(output2) <- NULL
colnames(output2) <- c("Sample","Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils","Basophils", "NeutrophilsPct", "LymphocytesPct", "MonocytesPct", "EosinophilsPct","BasophilsPct")
#Write to output entity
molgenis.addAll("${resultSetRepositoryName}", output2)
