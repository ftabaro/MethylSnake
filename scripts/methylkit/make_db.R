library(methylKit)

input_paths <- snakemake@input
methylkitdb_folder <- snakemake@params[["methylkitdb_folder"]]
methylkit_rdata_path <- snakemake@output[["methylRawObj"]]

sample_sheet <- read.csv(snakemake@input[["sample_sheet"]], header = TRUE)
sample.id <- sample_sheet$sample_name
treatment <- sample_sheet$treatment

message(sprintf("MethylKit database folder is: %s", methylkitdb_folder))
methylRawObj <- methRead(input_paths,
  sample.id=as.list(sample.id),
  assembly="hg19",
  dbtype="tabix",
  pipeline="bismarkCytosineReport",
  header=FALSE,
  skip=0,
  sep="\t",
  context="CpG",
  resolution="base",
  treatment=treatment,
  dbdir=methylkitdb_folder,
  mincov=10)
message("DB done.")

message(sprintf("MethylKit Rdata object path is: %s", methylkit_rdata_path))
saveRDS(methylRawObj, file = methylkit_rdata_path)
message("RData object saved to disk.")


pictures_folder <- snakemake@params[["pictures_folder"]]
for (i in seq_along(sample.id)) {
  sample <- sample.id[[i]]

  out.file <- file.path(pictures_folder, sprintf("%s_preliminary_stats.pdf", sample))
  message(sprintf("Sample id : %s - PDF path is: %s", sample, out.file))

  pdf(out.file, height=20, width=10)
  par(pty="s", mfrow=c(2,1))
  getMethylationStats(methylRawObj[[i]], plot=TRUE, both.strands=FALSE)
  getCoverageStats(methylRawObj[[i]], plot=TRUE, both.strands=FALSE)
  dev.off()

}

message("all done.")
