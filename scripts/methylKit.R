# R CMD BATCH '--args {alignments_dir} {methyldb_dir} {threads}' methylKit.R

library(methylKit)

args <- commandArgs(trailingOnly=TRUE)

alignments_dir <- args[1]
dbdir <- args[2]
ncores <- as.numeric(args[3])

file.list <- list.files(alignments_dir, pattern = "CpG.report.txt")
file.list <- as.list(file.list)

sample.id <- gsub("^(.+)_1_val_1.+$", "\\1", basename(file.list))
treatment <- ifelse(grepl("SRR", sample.id), 0, 1)

methylRawObj <- methRead(file.list, sample.id=sample.id, assembly="hg19", treatment=treatment, context="CpG", dbtype="tabix", dbdir=dbdir)

for (i in seq_along(methylRawObj)) {
  sample <- sample.id[i]
  pdf(sprintf("pictures/%s.pdf", sample), height=20, width=10)
  par(pty="s", mfrow=c(2,1))
  getMethylationStats(methylRawObj[[i]], plot=TRUE, both.strands=FALSE)
  getCoverageStats(methylRawObj[[i]], plot=TRUE, both.strands=FALSE)
  dev.off()
}

meth <- unite(methylRawObj)
pdf(sprintf("pictures/diagnostics.pdf", sample), height=10, width=10)
par(pty="s", mfrow=c(2,2))
getCorrelation(meth, plot = TRUE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth)
PCASamples(meth, screeplot=TRUE)
dev.off()

tiles <- tileMethylCounts(methylRawObj, win.size=500, step.size=500)
diffMeth <- calculateDiffMeth(tiles, overdispersion="MN", test="Chisq", mc.cores=ncores)

diff25p.hyper <- getMethyDiff(diffMeth, difference=25, qvalue=0.01, type="hyper")
diff25p.hypo <- getMethyDiff(diffMeth, difference=25, qvalue=0.01, type="hypo")

exclude <- diffMeth$chr[grepl("chr.*_", diffMeth$chr)]
pdf("pictures/differentiallyMethylatedPerChromosome.pdf")
diffMethPerChr(diffMeth,plot=TRUE, qvalue.cutoff=0.01, meth.cutoff=25, exclude = exclude)
dev.off()



