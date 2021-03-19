### load the library
library(ATACseqQC)
### input the bamFile from the ATACseqQC package
bamfile <- "~/mnt/bioinformaticsRTP/BlanksData/BulkATACseq/alignment/Control_1h_2.bam"
bamfile.labels <- gsub(".bam", "", basename(bamfile))
outFolder <- "atacseqqc"
dir.create(outFolder, recursive=TRUE)

### IGV snapshot function:
#source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))

png(paste0(outFolder, "/", "atacSeqQc_01_libraryComplexity_",bamfile.labels,".png"))
estimateLibComplexity(readsDupFreq(bamfile))
dev.off()

### generate fragment size distribution
png(paste0(outFolder, "/", "atacSeqQc_02_fragmentSizeDistribution_",bamfile.labels,".png"))
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()

## bamfile tags to be read in
possibleTag <- list(
  "integer"=c(
    "AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2",
    "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
    "TC", "UQ"),
  "character"=c(
    "BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
    "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
    "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
    "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
    "U2")
)
library(Rsamtools)
bamTop100 <- scanBam(
  BamFile(bamfile, yieldSize = 100),
  param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
#tags
## shift the coordinates of 5'ends of alignments in the bam file
library(BSgenome.Hsapiens.NCBI.GRCh38)
#seqlev <- "1" ## subsample data for quick run
seqlev <- c(1:22,"X","Y")
which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
outBamFl <- paste0(outFolder, "/", "atacSeqQc_03_bamshift_",bamfile.labels,".bam")
gal1 <- shiftGAlignmentsList(gal, outbam=outBamFl)

### PT score is calculated as the coverage of promoter divided by the coverage
### of its transcript body. PT score will show if the signal is enriched in
### promoters.
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevelsStyle(txdb) <- "NCBI"
txs <- transcripts(txdb)
pt <- PTscore(gal1, txs)
png(paste0(outFolder, "/", "atacSeqQc_04_fragmentSizeDistribution_",bamfile.labels,".png"))
plot(pt$log2meanCoverage, pt$PT_score,
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript",
     main=paste0("Transcript Vs Promoter Counts: ",bamfile.labels),
     col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.1)
)
abline(h=0, col=2, lty=2)
dev.off()


### Nucleosome Free Regions (NFR) score
nfr <- NFRscore(gal1, txs)
png(paste0(outFolder, "/", "atacSeqQc_05_NFRscore_",bamfile.labels,".png"))
plot(nfr$log2meanCoverage, nfr$NFR_score,
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main=paste0("NFRscore for 200bp flanking TSSs: ",bamfile.labels),
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()

### Transcription Start Site (TSS) Enrichment Score
tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore
png(paste0(outFolder, "/", "atacSeqQc_06_TSSscore_",bamfile.labels,".png"))
plot(100*(-9:10-.5), tsse$values, type="b",
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()

### Split reads:
library(phastCons100way.UCSC.hg38)
## run program for chromosome 1 only
txs <- txs[seqnames(txs) %in% "22"]
#txs <- txs[seqnames(txs) %in% c(1:21,"X","Y")]
genome <- Hsapiens
## split the reads into NucleosomeFree, mononucleosome,
## dinucleosome and trinucleosome.
## and save the binned alignments into bam files.
outPath <- paste0(outFolder, "/", "atacSeqQc_split_", bamfile.labels)
objs <- splitGAlignmentsByCut(
  gal1, txs=txs, genome=genome, outPath = outPath,
  conservation=phastCons100way.UCSC.hg38)
#dir(outPath)

### Heatmap and coverage curve for nucleosome positions
library(ChIPpeakAnno)
bamfiles <- file.path(
  outPath,
  c(
    "NucleosomeFree.bam",
    "mononucleosome.bam",
    "dinucleosome.bam",
    "trinucleosome.bam"
  )
)
## Plot the cumulative percentage of tag allocation in nucleosome-free
## and mononucleosome bam files.
png(paste0(outFolder, "/", "atacSeqQc_07_NucFree_",bamfile.labels,".png"))
cumulativePercentage(bamfiles[1:2], as(seqinfo(Hsapiens)["1"], "GRanges"))
dev.off()

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
## estimate the library size for normalization
librarySize <- estLibSize(bamfiles)
## calculate the signals around TSSs.
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(
  gal=objs[c("NucleosomeFree",
    "mononucleosome",
    "dinucleosome",
    "trinucleosome")],
  TSS=TSS,
  librarySize=librarySize,
  seqlev=seqlev,
  TSS.filter=0.5,
  n.tile = NTILE,
  upstream = ups,
  downstream = dws
)
## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
#plot heatmap
png(paste0(outFolder, "/", "atacSeqQc_08_NucleosomeHeatmap_",bamfile.labels,".png"))
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)
dev.off()

## get signals normalized for nucleosome-free and nucleosome-bound regions.
png(paste0(outFolder, "/", "atacSeqQc_09_NucleosomeOcc_",bamfile.labels,".png"))
out <- featureAlignedDistribution(
  sigs,
  reCenterPeaks(TSS, width=ups+dws),
  zeroAt=.5, n.tile=NTILE, type="l",
  ylab="Averaged coverage"
)
## rescale the nucleosome-free and nucleosome signals to 0~1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n",
        xlab="Position (bp)",
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1,
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()

library(MotifDb)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
# print(CTCF[[1]], digits=2)
png(paste0(outFolder, "/", "atacSeqQc_10_CTCFmotifLocalDp_",bamfile.labels,".png"))
sigs <- factorFootprints(
  outBamFl, pfm=CTCF[[1]],
  genome=genome,
  min.score="90%", seqlev=seqlev,
  upstream=100, downstream=100
)
dev.off()

### CTCF V-plot
png(paste0(outFolder, "/", "atacSeqQc_11_CTCFvplot1_",bamfile.labels,".png"))
vp <- vPlot(
  outBamFl, pfm=CTCF[[1]],
  genome=genome, min.score="90%", seqlev=seqlev,
  upstream=200, downstream=200,
  ylim=c(30, 250), bandwidth=c(2, 1)
)
dev.off()

png(paste0(outFolder, "/", "atacSeqQc_12_CTCFmotifvplot2_",bamfile.labels,".png"))
distanceDyad(vp, pch=20, cex=.5)
dev.off()

#####
#####
#####
