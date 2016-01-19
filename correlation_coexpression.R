##########################################################
#Author: Marije van der Geest                            #
#                                                        #
#Generates line plot for gluten specific T-cell          #
#expression per gene under different conditions          #
#(timpoints 0, 10, 30, 180)                              #
##########################################################

########
# Libs #
########
library(gplots)

#Make connection
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")

#Get and process genes (input)
genes <- "${genes}"
genes <- unlist(strsplit(genes, '[,]'))

#Create query for filtering data
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("GeneNames==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Get data
LLDeep <- molgenis.get("[INPUTDATA]", num=100000, q=q)
rownames(LLDeep) <- LLDeep[,1]
LLDeep <- LLDeep[,-1]
genes.length <- nrow(LLDeep)


#Get query for gene id -> gene name conversion
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Convert genes
conversion.table <- molgenis.get("GeneInfo", q=q)
row.names(LLDeep) = as.character(conversion.table[do.call(c,lapply(row.names(LLDeep), function(x) which(as.character(conversion.table[,1]) == x))),2])

#Perform correlation (Spearman)
correlation.matrix <- cor(t(LLDeep), method="spearman")

#Create image
png("${outputFile}")
colors <- colorRampPalette(c("gold", "midnightblue", 'blue'))(10)
heatmap.2(correlation.matrix, trace="none", 
          colsep=c(1:genes.length), rowsep=c(1:genes.length), 
          sepcolor="snow2", sepwidth=c(0.05,0.05), col=colors, 
          cexRow=1.2, cexCol=1.2, key.par=list(c(-1, 0, 0.5, 1)), margin=c(7,10))