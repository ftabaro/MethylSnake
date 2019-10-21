library(methylKit)
library(genomation)
library(GenomicRanges)

dmr <- readRDS(snakemake@input[["dmr"]])
annotation_file <- snakemake@input[["annotation"]]

gene.obj <- readTranscriptFeatures(annotation_file)
seqlevelsStyle(gene.obj) <- "UCSC"

dmr_gr <- as(dmr, "GRanges")

seqlevels(dmr_gr, pruning.mode="coarse") <- seqlevels(gene.obj)

annotated <- annotateWithGeneParts(dmr_gr, gene.obj)
saveRDS(annotated, snakemake@output[["dmr_annotated"]])

p <- getTargetAnnotationStats(annotated, percentage=TRUE, precedence=TRUE)
cat(p, file = snakemake@output[["genomic_features_percentages"]])


pdf(snakemake@output[["gene_part_annotation_plot"]])
plotTargetAnnotation(annotated, precedence = TRUE, main="Annotation of DMRs")
dev.off()

tbl <- getAssociationWithTSS(annotated)
write.csv(tbl, file=snakemake@output[["tss_association"]])

