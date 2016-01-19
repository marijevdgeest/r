##########################################################
#Authors: Marije van der Geest                           #
#                                                        #
#Generates heatmap for gamma/delta T-cell expression     #
#per gene under different conditio                       #
##########################################################

########
# Libs #
########
library("ggplot2")
library("RColorBrewer")
require("reshape2")


#Create connection with molgenis
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")

#Get and process genes (input)
id.my <- "${genes}"
genes <- unlist(strsplit(id.my, '[,]'))

#Get sample info
samples <- molgenis.get("[SAMPLEDATA]")
rownames(samples) <- samples[,1]
samples <- samples[,-1]

#Create query for filtering data
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("GenesIds==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Access data
gdTcells.means <- molgenis.get("[INPUTDATA]", q=q)
rownames(gdTcells.means) <- gdTcells.means[,1]
gdTcells.means <- gdTcells.means[,-1]

#Create query for gene id >- gene name conversion
qs <- ""
for (i in 1:length(genes)){
  qs[i] <- paste("EnsemblGeneID==", genes[i], sep="")
}
q=paste(qs, collapse=",")

#Convert genes
conversion.table <- molgenis.get("GeneInfo", q=q)
row.names(gdTcells.means) = as.character(conversion.table[do.call(c,lapply(row.names(gdTcells.means), function(x) which(as.character(conversion.table[,1]) == x))),2])

#Format data
gdTcells.means.df  <- melt(as.matrix(gdTcells.means))
#Combine vaues and sample info
betterSamples <- cbind(rownames(samples), samples)
colnames(betterSamples)[1] <- "Sname"
finaldf <- merge(gdTcells.means.df, betterSamples, by.x="Var2", by.y="Sname")
#Filter data for relevant conditions
finaldf <- finaldf[finaldf$State =="Active" | finaldf$State =="Control",]
finaldf <- finaldf[finaldf$Tissue =="IEL",]

#Set colors and breaks
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
b <- c(-1.5,0,1.5)

png("${outputFile}", width=650, height=500)

#Create heatmap
p <- ggplot(finaldf, aes(y=Var1, x=Stimulation, group=State)) + geom_tile(aes(fill = log10(value)), colour = "White") + geom_point(data = finaldf, aes(size="0", shape = NA), colour = "grey50")
p + facet_grid(CellType~State) + scale_fill_gradientn(limits = c(-3,3),colours=hmcol,breaks=b,labels=format(b)) + xlab("Stimulation") + ylab("Genes") + guides(size=guide_legend("-Inf", override.aes=list(shape=15, size = 10))) + theme_grey(base_size = 16) 