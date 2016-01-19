##########################################################
#Authors: Marije van der Geest                           #
#                                                        #
#Generates heatmap for rpkm expression data for          #
#different cell types                                    #
##########################################################

########
# Libs #
########
library("gplots")
library("fields")


#Make connection with Molgenis
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")

#Get and process genes (input)
genes <- strsplit("${genes}", '[,]')


#Query to filter data
query <- lapply(genes, FUN=function(geneName) paste("GeneNames==",geneName, collapse=" or ", sep=""))
dane <- molgenis.get("[INPUTDATA]", q=query)
#Format data
rownames(dane) <- dane[,1] 
dane <- dane[,2:8]
mysep <- c(1:nrow(dane))              


#Create image
colors <- c("#FFFFFF", "#EFCC8C", "#D7A069", "#C1755B", "#A63A4A")
par(oma=c(20,0,0,20))
png("${outputFile}", width=500, height=500)
heatmap.2(as.matrix(dane),col=colors, trace="none",margin=c(10,10),scale="none",dendrogram="none", Rowv="NULL", 
          Colv="NULL",key=F,breaks=c(-2,-0.1,1,5,10,max(dane)), cexRow=1.2,
          cexCol=1.2, lhei=c(1,8), rowsep=mysep, colsep=c(1:28), sepcolor="snow2", sepwidth=c(0.01,0.01),
          labCol=c("NK cells","B-cells","Monocytes","T-memory cells","CD4 T-cells","CD8 T-cells","Granulocytes"))

image.plot(matrix(-3:12,16,8),legend.only=TRUE, horizontal=TRUE, col=colors, legend.shrink=0.25, legend.width=0.5,
           axis.args=list(at=c(-3,0,3,6,9,12),labels=c("NA",0,3,5,10,"MAX")),smallplot=c(0.51,0.9,0.96,0.99))