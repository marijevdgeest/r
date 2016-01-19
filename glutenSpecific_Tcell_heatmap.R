##########################################################
#Authors: Yang Li & Marije van der Geest                 #
#                                                        #
#Generates heatmap for gluten specific T-cell expression #
#per gene under different conditions (timpoints 0, 10,   #
#30, 180)                                                #
##########################################################

########
# Libs #
########
library("gplots")
library("RColorBrewer")


#Create connection with molgenis
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")

#Get and process genes (input)
id.my <- "${genes}"
genes <- unlist(strsplit(id.my, '[,]'))

#Create query for filtering data
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("GeneIds==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Get data
groupMean <- molgenis.get("[INPUTDATA]", q=q)

#Format data
rownames(groupMean) <- groupMean[,1]
groupMean <- groupMean[,-1]
data.plot <- groupMean[genes,]

#Remove duplicates
genes.remove <- which(row.names(data.plot) != genes)
if (length(genes.remove) != 0){
  genes <- genes[-genes.remove]
  data.plot <- na.omit(data.plot)
}

#Create query for gene id -> gene name conversion
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")
#Convert genes
conversion.table <- molgenis.get("GeneInfo", q=q)
row.names(data.plot) = as.character(conversion.table[do.call(c,lapply(row.names(data.plot), function(x) which(as.character(conversion.table[,1]) == x))),2])


#Set colors
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#Create image
png("${outputFile}",width=500, height=400)
heatmap.2(as.matrix(data.plot), col = hmcol, trace="none", margin=c(5, 10),
          scale="row",dendrogram="row",Colv=FALSE,colsep=1:3,srtCol=0,
          sepcolor="gray", cexRow=1.2, cexCol=1.3, labCol=c("0","10","30","180"))