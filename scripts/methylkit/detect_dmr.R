library(methylKit)
library(GenomicRanges)

methylMergedObj <- readRDS(snakemake@input[["methylMergedObj"]])

suffix <- sprintf("tiled_ws%d_ss%d", snakemake@params[["window_size"]], snakemake@params[["step_size"]])
tiles <- tileMethylCounts(methylMergedObj,
  win.size=snakemake@params[["window_size"]],
  step.size=snakemake@params[["step_size"]],
  suffix=suffix,
  mc.cores=snakemake@threads)
saveRDS(tiles, snakemake@output[["tileCounts"]])

diffMeth <- calculateDiffMeth(tiles,
  overdispersion="MN", test="Chisq",
  mc.cores=snakemake@threads)
saveRDS(diffMeth, snakemake@output[["dmr"]])

diff.hyper <- getMethylDiff(diffMeth,
  difference=snakemake@params[["difference"]],
  qvalue=snakemake@params[["qvalue"]],
  type="hyper")
saveRDS(diff.hyper, snakemake@output[["dmrHyper"]])

diff.hypo <- getMethylDiff(diffMeth,
  difference=snakemake@params[["difference"]],
  qvalue=snakemake@params[["qvalue"]],
  type="hypo")
saveRDS(diff.hypo, snakemake@output[["dmrHypo"]])

dmrPerChr <- diffMethPerChr(diffMeth, plot = FALSE,
  qvalue.cutoff=snakemake@params[["qvalue"]],
  meth.cutoff=snakemake@params[["difference"]])
saveRDS(dmrPerChr, snakemake@output[["dmrPerChromosome"]])

chroms <- levels(seqnames(as(diffMeth, "GRanges")))
exclude <- as.factor(grep("_", chroms, value = TRUE))

pdf(snakemake@output[["dmrPerChromosome_plot"]])
diffMethPerChr(diffMeth,plot = TRUE,
  qvalue.cutoff=snakemake@params[["qvalue"]],
  meth.cutoff=snakemake@params[["difference"]],
  exclude = exclude)
dev.off()



