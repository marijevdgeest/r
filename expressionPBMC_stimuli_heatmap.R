##########################################################
#Author: Marije van der Geest                            #
#                                                        #
#Generates heatmap for the gene expression in PBMC       #
#samples for different stimuli experiments.              #
##########################################################

########
# Libs #
########
library("fields")
library("gplots")


#Create connection to molgenis
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")

#Get and process genes (input)
gene_names <- "${genes}"
genes <- unlist(strsplit(gene_names, '[,]'))


#Create query for filtering data
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("GeneIds==", genes[i], sep="")
}
q=paste(qs, collapse=",")
#Get data
rpkm.mean.mat <- molgenis.get("[INPUTDATA]", q=q)


#Format data.
rownames(rpkm.mean.mat) <- rpkm.mean.mat[,1]
rpkm.mean.mat <- rpkm.mean.mat[,-1]
dane <- rpkm.mean.mat[genes,]

genes.remove <- which(row.names(dane) != genes)
if (length(genes.remove) != 0){
  genes <- genes[-genes.remove]
  dane<- na.omit(dane)
}

#Create query for gene id -> gene name conversion. 
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")
#Convert genes
conversion.table <- molgenis.get("GeneInfo", q=q)
row.names(dane) = as.character(conversion.table[do.call(c,lapply(row.names(dane), function(x) which(as.character(conversion.table[,1]) == x))),2])


#Order stimuli
dane <- dane[,order(names(dane), decreasing=T)]


#create image
colors <- c("#FFFFFF", "#EFCC8C", "#D7A069", "#C1755B", "#A63A4A")
png("${outputFile}",width=500, height=500)
heatmap.2(as.matrix(dane), trace="none", col=colors, dendrogram="none", margin=c(15,10), cexCol=1.2, cexRow=1.2, breaks=c(-2,-0.1,1,5,10,max(dane, na.rm=T)),
          rowsep=c(1:length(genes)), colsep=c(1:15), sepcolor="snow2", sepwidth=c(0.01,0.01), key=F, lhei=c(1,4), Rowv="NULL", Colv="NULL", srtCol=45)
image.plot(matrix(-3:12,16,1),legend.only=TRUE, horizontal=TRUE, col=colors, legend.shrink=0.25, legend.width=0.5,
           axis.args=list(at=c(-3,0,3,6,9,12),labels=c("NA",0,3,5,10,"MAX")),smallplot=c(0.51,0.9,0.96,0.99))

