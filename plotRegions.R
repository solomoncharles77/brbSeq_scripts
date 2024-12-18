
library(karyoploteR)

TP53.region <- toGRanges("chr17:7,564,422-7,602,719")

kp <- plotKaryotype(zoom = TP53.region)


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)

kp <- plotKaryotype(zoom = TP53.region)
kpPlotGenes(kp, data=genes.data)


kp <- plotKaryotype(zoom = TP53.region, cex=2)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)



testRegion <- toGRanges("chr9:40,457,292-40,468,587")
testRegion <- toGRanges("chr9:39e6-41e6")
testRegion <- toGRanges("chr9:100,887,633-100,887,633")
tp <- plotKaryotype(zoom = testRegion)

testRegion <- toGRanges("chr9:39e6-40e6")
tp <- plotKaryotype(genome="hg38", zoom = testRegion)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=tp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(tp, data=genes.data)



library(TxDb.Hsapiens.UCSC.hg19.knownGene)

kp <- plotKaryotype(plot.type=4, zoom="chr1:68e6-73e6")
kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e6)
kpAddLabels(kp, labels = "Trait 1", srt=90, pos=3, cex=1.8, label.margin = 0.025)
kpAxis(kp, ymin=0, ymax=10, r0=0.2)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set3", ymax=10, r0=0.2)

kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3, cex=1.6, col="red", r0=0.2)
kpPoints(kp, data = top.snps, pch=1, cex=1.6, col="red", lwd=2, ymax=10, r0=0.2)

genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8, gene.name.position = "left")
