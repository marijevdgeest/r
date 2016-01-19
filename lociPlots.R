##########################################################
#Authors: Patrick Deelen & Marije van der Geest          #
#                                                        #
#Loci plots    (STILL IN DEVELOPMENT)                    #
##########################################################

#############
# Resources #
#############
# chrLoci <- "chr19"
# startLoci <- 48956172
# endLoci <- 49456172
# gsnp <- "rs516246"

#Get input data
chrLoci <- "${chrLoci}"
startLoci <- "${startLoci}"
endLoci <- "${endLoci}"
gsnp <- "${gsnp}"


########
# Libs #
########
library(gplots)

#Make connection
source("http://localhost:8080/molgenis.R")
molgenis.login("[USERNAME]", "[PASSWORD]")


#Create query for filtering data
q1 <- paste("chrom==", paste("chr",chrLoci,sep=""), sep="")

lociTranscripts <- molgenis.get("[INPUTDATA]", q=q1)

transcripts <-  lociTranscripts[lociTranscripts$txStart>=startLoci & lociTranscripts$txStart<=endLoci,]
rownames(transcripts) <- transcripts[,1]
transcripts <- transcripts[,-1]

q2 <- paste("gSNP==", gsnp, sep="")
locusEqtlGenes <- molgenis.get("[INPUTDATA]", q=q2)
q3 <- paste("Chromosome==", chrLoci, sep="")
snpsToPlot <- molgenis.get("[INPUTDATA]t", q=q3)
sourceTranslate <- molgenis.get("[INPUTDATA]")
sourceTranslate <- sourceTranslate[,-1]
rownames(sourceTranslate) <- sourceTranslate[,1]
sourceTranslate <- as.data.frame(sourceTranslate[,-1], row.names=as.character(sourceTranslate[,1]), stringsAsFactors = FALSE)


locusSnpsToPlot <- snpsToPlot [snpsToPlot$POS >= startLoci & snpsToPlot$POS <= endLoci, ]
transcripts$source2 <- sourceTranslate[as.character(transcripts$source),]


geneStarts <- aggregate(txStart ~ name2, data = transcripts, FUN = min)
geneStops <- aggregate(txEnd ~ name2, data = transcripts, FUN = max)

selectedTranscriptsNames <- character(length = length(levels(transcripts$name2)))
i <- 1;

for(gene in levels(transcripts$name2)){
  geneTranscipts <- transcripts[transcripts$name2 == gene,]
  
  a <- match("protein coding", geneTranscipts$source2 )
  if(is.na(a)){
    a <- 1
  }
  
  selectedTranscriptsNames[i] <- as.character(geneTranscipts$name)[a]
  i <- i + 1
  
}


genes <- transcripts[selectedTranscriptsNames %in% transcripts$name, c("name2", "value", "source", "source2", "strand")]
genes <- merge(genes, geneStarts)
genes <- merge(genes, geneStops)

genes <- genes[order(genes$txStart),]

geneCol <- c("protein coding" = "dodgerblue", "long non-coding" = "red3", "other" = "grey50", "small non-coding" = "forestgreen", "pseudogene" = "orange")

png("${outputFile}", width=650, height=650)
par(xpd=NA, mar = c(3,0,0.5,1))
plot.new()
plot.window(xlim = c(as.numeric(startLoci), as.numeric(endLoci)), ylim = c(0,3))

a <- TRUE
b <- TRUE


for(i in 1:nrow(genes)){
  gene <- genes[i,]
  
  y1 <- 0.5
  y2 <- 0.7
  if(gene$strand == "-"){
    y1 <- y1 + 1
    y2 <- y2 + 1
    if (a){
      a <- FALSE
      y1 <- y1 + 0.5
      y2 <- y2 + 0.5
    } else {
      a <- TRUE
    }
  } else if (b){
    b <- FALSE
    y1 <- y1 + 0.5
    y2 <- y2 + 0.5
  } else {
    b <- TRUE
  }
  
  
  geneStart <- gene$txStart
  geneEnd <- gene$txEnd
  
  rect(geneStart, y1, geneEnd, y2, col = geneCol[as.character(gene$source2)])
  
  if(gene$strand == "-"){
    lines(x = c(geneEnd, geneEnd), y = c(y2, y2 + 0.1))
    arrows(geneEnd, y2 + 0.1, geneStart, y2 + 0.1, length = 0.05)
  } else {
    lines(x = c(geneStart, geneStart), y = c(y2, y2 + 0.1))
    arrows(geneStart, y2 + 0.1, geneEnd, y2 + 0.1, length = 0.05)
  }
  
  if(gene$name2 %in% locusEqtlGenes$Probe)
  {
    if(gene$strand == "-"){
      text(x = geneStart + ((geneEnd - geneStart)/2), y = 0.5 + 1.9, labels = gene$value, srt = 45, adj = c(0,0.5), cex = 0.7)
    } else {
      text(x = geneStart + ((geneEnd - geneStart)/2), y = 0.5 - 0.1, labels = gene$value, srt = 45, adj = c(1,0.5), cex = 0.7)
      
    }
  }
  
}

snpPch <- rep(3, times = nrow(locusSnpsToPlot))
snpPch[ as.character(locusSnpsToPlot$Snp_Top_Proxy) == "TopSNP" ] <- 7
snpPch[ as.character(locusSnpsToPlot$Snp_Top_Proxy) == "TopSNP_noeQTL" ] <- 4


snpCex <- rep(0.5, times = nrow(locusSnpsToPlot))
snpCex[ as.character(locusSnpsToPlot$Snp_Top_Proxy) == "TopSNP" ] <- 1.5
snpCex[ as.character(locusSnpsToPlot$Snp_Top_Proxy) == "TopSNP_noeQTL" ] <- 1.5

points(locusSnpsToPlot$POS, y = rep(1.4, times = nrow(locusSnpsToPlot)), pch = snpPch, cex = snpCex, col = geneCol)




legend(x = "bottom", legend = names(geneCol), fill = geneCol, ncol = 5, bty = "n", cex = 1)
legend(x = "bottom", legend = c("Top SNP", "Top SNP, no eQTL", "Proxy SNP"), pch = c(7,4,3), ncol = 5, bty = "n", cex = 1, inset = c(0,-0.05))


