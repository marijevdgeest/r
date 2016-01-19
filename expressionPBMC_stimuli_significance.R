##########################################################
#Author: Marije van der Geest                            #
#                                                        #
#Generates heatmap for significance of stimulated PMBC   #
#expression.                                             #
##########################################################

########
# Libs #
########
library("gplots")


#Create connection to molgenis
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")

#Get and process genes (input)
genes.id <- "${genes}"
genes <- unlist(strsplit(genes.id, '[,]'))

#Create query for filtering data
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("gene_id==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Get data
stimuli.significance <- molgenis.get("[INPUTDATA]", q=q)
rownames(stimuli.significance) <- stimuli.significance[,1]
#Format data
genes.selected.list <- lapply(genes, function(x) data.frame(stimuli.significance[x,]))
genes.selected <- do.call(rbind, genes.selected.list)
genes.selected.plot <- na.omit(genes.selected[,4:15])


#Add missing genes with '0' (not significant under any condition)
for (gene in genes){
  if (!(gene %in% row.names(genes.selected.plot))) {genes.selected.plot[gene,]=0}}


#Create query for gene id -> gene name conversion
qs <- ""
for (i in 1:length(row.names(genes.selected.plot))){
  qs[i] <- paste("EnsemblGeneID==", row.names(genes.selected.plot)[i], sep="")
}
q=paste(qs, collapse=",")
#Convert genes
conversion.table <- molgenis.get("GeneInfo", q=q)
row.names(conversion.table) <- conversion.table[,1]

#Remove duplicate genes
for (x in row.names(genes.selected.plot)){
  if (!(x %in% row.names(conversion.table))){
    drop <- which(row.names(genes.selected.plot) == x)
    genes.selected.plot <- genes.selected.plot[-drop,]
  }
}

row.names(genes.selected.plot) = as.character(conversion.table[do.call(c,lapply(row.names(genes.selected.plot), function(x) which(as.character(conversion.table[,1]) == x))),2])


#Order stimuli conditions
genes.selected.plot.ordered <- genes.selected.plot[,order(names(genes.selected.plot), decreasing=T)]


#Create image
colors <- c("royalblue4", "red3")
png("${outputFile}",width=500, height=500)
heatmap.2(as.matrix(genes.selected.plot.ordered), cexCol=1.2, cexRow=1.2, col=colors, scale="none", 
          trace="none", key=F, dendrogram="none", srtCol=45, Colv="NULL", Rowv="NULL", margin=c(15,10), 
          sepcolor="snow2", sepwidth=c(0.01,0.01), rowsep=c(1:length(genes)), colsep=c(1:11))




#Generate legend
legend(x="topright", legend=c("Not significant", "Significant"), col=colors, pch=15, cex=1.2)