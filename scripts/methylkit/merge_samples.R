library(methylKit)

methylRawObj <- readRDS(snakemake@input[["methylRawObj"]])

meth <- unite(methylRawObj)
saveRDS(meth, snakemake@output[["methylMergedObj"]])

pdf(snakemake@output[["correlogram"]], height=30, width=30)
getCorrelation(meth, plot = TRUE)
dev.off()

pdf(snakemake@output[["clustering"]], height=7, width=21)
par(pty="s", mfrow=c(1,3))
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth)
PCASamples(meth, screeplot=TRUE)
dev.off()
