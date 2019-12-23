options(warn=1)
sink(stdout(), type = "message")
message("#################################\n## Starting methylKit analysis ##\n#################################")

suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))



######################
## READ CONFIG FILE ##
######################

# shared folders
methylkitdb_folder <- snakemake@config[["methylkitdb_folder"]]
alignments_folder  <- snakemake@config[["alignments_folder"]]

rdata_folder    <- snakemake@config[["rdata_folder"]]
bed_folder      <- snakemake@config[["bed_folder"]]
pictures_folder <- snakemake@config[["pictures_folder"]]
tables_folder   <- snakemake@config[["tables_folder"]]
message(sprintf("Folder summary:\n\tmethylkitdb_folder: %s\n\talignments_folder: %s\n\trdata_folder: %s\n\tpictures_folder: %s\n\ttables_folder: %s",methylkitdb_folder, alignments_folder, rdata_folder, pictures_folder, tables_folder))

# dmr paramters
window_size   <- as.numeric(snakemake@config[["dmr_window_size"]])
step_size     <- as.numeric(snakemake@config[["dmr_step_size"]])
diff          <- as.numeric(snakemake@config[["dmr_difference"]])
qvalue        <- as.numeric(snakemake@config[["dmr_qvalue"]])
min_per_group <- as.numeric(snakemake@config[["min_per_group"]])
message(sprintf("DMR calling parameters: ws=%d (%s) - ss=%d  (%s) - diff = %d (%s) - q-value=%.6f (%s) - min-per-group=%d (%s)",
  window_size, class(window_size),
  step_size, class(step_size),
  diff, class(diff),
  qvalue, class(qvalue),
  min_per_group, class(min_per_group)))

# annotation file
annotation_file <- snakemake@input[["annotation_file"]]

# other params
threads            <- snakemake@threads
genome_version     <- snakemake@config[["genome_version"]]
mate1_pattern      <- snakemake@config[["mate1_pattern"]]
samples_sheet_path <- snakemake@input[["sample_sheet"]]



#################
## BUILD PATHS ##
#################

# DMR and DMC paths
dmr_rdata_folder    <- file.path(rdata_folder, "dmr")
dmr_pictures_folder <- file.path(pictures_folder, "dmr")
dmr_tables_folder   <- file.path(tables_folder, "dmr")

dmc_rdata_folder    <- file.path(rdata_folder, "dmc")
dmc_pictures_folder <- file.path(pictures_folder, "dmc")
dmc_tables_folder   <- file.path(tables_folder, "dmc")

for(d in c(dmr_rdata_folder, dmr_pictures_folder, dmr_tables_folder, dmc_rdata_folder, dmc_pictures_folder, dmc_tables_folder)) {
  exit_code <- dir.create(d, showWarnings=FALSE, recursive=TRUE)
  if (exit_code) {
    message(sprintf("Created %s", d))
  }
}

# RData objects paths
methylRawObj_rds_path    <- file.path(rdata_folder, "methylRawObj.rds")
methylMergedObj_rds_path <- file.path(rdata_folder, "methylMergedObj.rds")
tileCounts_rds_path      <- file.path(rdata_folder, "tilesCounts.rds")
diffMeth_rds_path        <- file.path(rdata_folder, "methylDiffObj.rds")
diffMethTiles_rds_path   <- file.path(rdata_folder, "methylDiffTilesObj.rds")

dmr_paths <- c(
  rds_path                        = file.path(dmr_rdata_folder, "dmr.rds"),
  annotated_rds_path              = file.path(dmr_rdata_folder, "dmr_annotated.rds"),
  perChr_rds_path                 = file.path(dmr_rdata_folder, "dmrPerChromosome.rds"),
  gene_part_annotation_plot       = file.path(dmr_pictures_folder, "annotation.pdf"),
  perChr_plot_path                = file.path(dmr_pictures_folder, "dmrPerChromosome.pdf"),
  tss_association_table           = file.path(dmr_tables_folder, "dmr_tss_association.csv"),

  hyper_rds_path                  = file.path(dmr_rdata_folder, "hyper.rds"),
  hyper_annotated_rds_path        = file.path(dmr_rdata_folder, "hyper_dmr_annotated.rds"),
  hyperPerChr_rds_path            = file.path(dmr_rdata_folder, "hyperDmrPerChromosome.rds"),
  hyper_gene_part_annotation_plot = file.path(dmr_pictures_folder, "hyper_annotation.pdf"),
  hyperPerChr_plot_path           = file.path(dmr_pictures_folder, "hyperDmrPerChromosome.pdf"),
  hyper_tss_association_table     = file.path(dmr_tables_folder, "hyperDmr_tss_association.csv"),

  hypo_rds_path                   = file.path(dmr_rdata_folder, "hypo.rds"),
  hypo_annotated_rds_path         = file.path(dmr_rdata_folder, "hypo_dmr_annotated.rds"),
  hypoPerChr_rds_path             = file.path(dmr_rdata_folder, "hypoDmrPerChromosome.rds"),
  hypo_gene_part_annotation_plot  = file.path(dmr_pictures_folder, "hypo_annotation.pdf"),
  hypoPerChr_plot_path            = file.path(dmr_pictures_folder, "hypoDmrPerChromosome.pdf"),
  hypo_tss_association_table      = file.path(dmr_tables_folder, "hypoDmr_tss_association.csv")

)

dmc_paths <- c(
  rds_path                        = file.path(dmc_rdata_folder, "dmc.rds"),
  annotated_rds_path              = file.path(dmc_rdata_folder, "dmc_annotated.rds"),
  perChr_rds_path                 = file.path(dmc_rdata_folder, "dmcPerChromosome.rds"),
  gene_part_annotation_plot       = file.path(dmc_pictures_folder, "annotation.pdf"),
  perChr_plot_path                = file.path(dmc_pictures_folder, "dmcPerChromosome.pdf"),
  tss_association_table           = file.path(dmc_tables_folder, "dmc_tss_association.csv"),

  hyper_rds_path                  = file.path(dmc_rdata_folder, "hyper.rds"),
  hyper_annotated_rds_path        = file.path(dmc_rdata_folder, "hyper_dmc_annotated.rds"),
  hyperPerChr_rds_path            = file.path(dmc_rdata_folder, "hyperDmcPerChromosome.rds"),
  hyper_gene_part_annotation_plot = file.path(dmc_pictures_folder, "hyper_annotation.pdf"),
  hyperPerChr_plot_path           = file.path(dmc_pictures_folder, "hyperDmcPerChromosome.pdf"),
  hyper_tss_association_table     = file.path(dmc_tables_folder, "hyperDmc_tss_association.csv"),

  hypo_rds_path                   = file.path(dmc_rdata_folder, "hypo.rds"),
  hypo_annotated_rds_path         = file.path(dmc_rdata_folder, "hypo_dmc_annotated.rds"),
  hypoPerChr_rds_path             = file.path(dmc_rdata_folder, "hypoDmcPerChromosome.rds"),
  hypo_gene_part_annotation_plot  = file.path(dmc_pictures_folder, "hypo_annotation.pdf"),
  hypoPerChr_plot_path            = file.path(dmc_pictures_folder, "hypoDmcPerChromosome.pdf"),
  hypo_tss_association_table      = file.path(dmc_tables_folder, "hypoDmc_tss_association.csv")

)

# plot paths
correlogram_path <- file.path(pictures_folder, "samples_correlation.pdf")
clustering_path  <- file.path(pictures_folder, "samples_clustering_pca.pdf")



##################################
## INIT ANNOTATION GRANGES LIST ##
##################################

# initialize annotation object
message(sprintf("Using annotation file from %s", annotation_file))
gene.obj                 <- genomation::readTranscriptFeatures(annotation_file)
seqlevelsStyle(gene.obj) <- "UCSC"



######################
## HELPER FUNCTIONS ##
######################

make_methylkitdb_object <- function (input_paths, sample.id, treatment, methylRawObj_rds_path, .methylkitdb_folder = methylkitdb_folder) {

  message(sprintf("Input paths (%d - %s): ", length(input_paths), class(input_paths)), paste(input_paths, collapse=", "))
  message(sprintf("Samples ids (%d -  %s): ", length(sample.id), class(sample.id)), paste(sample.id, collapse=", "))
  message(sprintf("Treatment vector (%d - %s): ", length(treatment), class(treatment)), paste(treatment, collapse=", "))

  if (is.null(.methylkitdb_folder)) {
    # work in memory
    methylRawObj <- methRead(as.list(input_paths),
      sample.id=as.list(sample.id),
      assembly=genome_version,
      dbtype=NA,
      pipeline="bismarkCytosineReport",
      header=FALSE,
      skip=0,
      sep="\t",
      context="CpG",
      resolution="base",
      treatment=treatment,
      mincov=10)
  }else{
    message("MethylKit database folder is: ", .methylkitdb_folder)

    # write tabix files to disk (SLOW)
    methylRawObj <- methRead(as.list(input_paths),
      sample.id=as.list(sample.id),
      assembly=genome_version,
      dbtype="tabix",
      pipeline="bismarkCytosineReport",
      header=FALSE,
      skip=0,
      sep="\t",
      context="CpG",
      resolution="base",
      treatment=treatment,
      dbdir=.methylkitdb_folder,
      mincov=10)
    message("DB done.")
  }

  message(sprintf("MethylKit Rdata object path is: %s", methylRawObj_rds_path))
  saveRDS(methylRawObj, file = methylRawObj_rds_path)
  message("RData object saved to disk.")

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

  message("MethylRawObj successfully computed.")
  return(methylRawObj)
}



merge_methylRawObj <- function(methylRawObj, .methylMergedObj_rds_path = methylMergedObj_rds_path, .correlogram_path = correlogram_path, .clustering_path = clustering_path) {

    meth <- unite(methylRawObj, mc.cores=threads, min.per.group=as.integer(min_per_group))
    saveRDS(meth, .methylMergedObj_rds_path)

    pdf(.correlogram_path, height=30, width=30)
    getCorrelation(meth, method="pearson", plot = TRUE)
    dev.off()

    pdf(.clustering_path, height=7, width=21)
    par(pty="s", mfrow=c(1,3))
    clusterSamples(meth,
      dist="correlation",
      method="ward",
      sd.threshold=0.5,
      filterByQuantile=TRUE,
      sd.filter=TRUE,
      plot=TRUE)
    PCASamples(meth)
    PCASamples(meth, screeplot=TRUE)
    dev.off()

    message("Merging complete.")
    return(meth)
}



compute_tiles <- function(methylMergedObj, ws = window_size, ss = step_size, .tileCounts_rds_path = tileCounts_rds_path) {

  message(sprintf("Computing tiles: ws=%d - ss=%d", ws, ss))

  suffix <- sprintf("tiled_ws%d_ss%d", ws, ss)

  tiles <- tileMethylCounts(methylMergedObj,
    win.size=ws,
    step.size=ss,
    suffix=suffix,
    mc.cores=threads)

  saveRDS(tiles, .tileCounts_rds_path)
  message("Tiling completed.")
  return(tiles)
}



add_number_to_path <- function(path, i) {
  sub("(.*)\\.(.*)$", sprintf("\\1\\.%d.\\2", i), path)
}



diff_meth_analysis <- function (methDiffObj,
  dmr_path, dmrPerChr_path, dmrPerChr_plot,
  dmr_annotated_path, gene_part_annot_plot, tss_association,
  direction = "all") {

  if (is.null(direction)){
    direction <- "all"
  }

  dmr <- getMethylDiff(methDiffObj,
    difference=diff,
    qvalue=qvalue,
    type=direction)
  saveRDS(dmr, dmr_path)

  if (direction == "all"){
    dmrPerChr <- diffMethPerChr(methDiffObj,
      plot = FALSE,
      qvalue.cutoff=qvalue,
      meth.cutoff=diff )
    saveRDS(dmrPerChr, dmrPerChr_path)

    chroms <- levels(seqnames(as(methDiffObj, "GRanges")))
    exclude <- as.factor(grep("_", chroms, value = TRUE))

    pdf(dmrPerChr_plot)
    diffMethPerChr(methDiffObj, plot = TRUE,
      qvalue.cutoff=qvalue,
      meth.cutoff=diff,
      exclude = exclude)
    dev.off()
  }

  message("Differential methylation analysis complete.")

  if(nrow(dmr) > 0){
    message(sprintf("Detected %d DMR (type = %s).", nrow(dmr), direction))

    dmr_gr <- as(dmr, "GRanges")

    seqlevels(dmr_gr, pruning.mode="coarse") <- seqlevels(gene.obj)

    annotated <- annotateWithGeneParts(dmr_gr, gene.obj)
    saveRDS(annotated, dmr_annotated_path)

    pdf(gene_part_annot_plot)
    plotTargetAnnotation(annotated,
      precedence = TRUE,
      main="Annotation of DMRs")
    dev.off()

    tbl <- getAssociationWithTSS(annotated)
    write.csv(tbl, file=tss_association)

    # cat(tbl, "\n")
    message("Annotation complete.")
  } else {
    warning("NO DMR TO ANNOTATE")
    saveRDS("NO DMR TO ANNOTATE", dmr_annotated_path)
  }

}

compute_diffMeth <- function(methylMergedObj, diffMethObj_rds_path = diffMeth_rds_path) {

    diffMeth <- calculateDiffMeth(methylMergedObj,
      overdispersion="MN",
      test="Chisq",
      mc.cores=threads)
    saveRDS(diffMeth, diffMethObj_rds_path)

    message("DiffMeth object computed.")
    return(diffMeth)
}

do_diff_meth_analysis <- function(diffMethObj, path_vector) {

  diff_meth_analysis(diffMethObj,
    dmr_path = path_vector["rds_path"],
    dmrPerChr_path = path_vector["perChr_rds_path"],
    dmrPerChr_plot = path_vector["perChr_plot_path"],
    dmr_annotated_path = path_vector["annotated_rds_path"],
    gene_part_annot_plot = path_vector["gene_part_annotation_plot"],
    tss_association = path_vector["tss_association_table"])

  diff_meth_analysis(diffMethObj, direction="hypo",
    dmr_path = path_vector["hypo_rds_path"],
    dmrPerChr_path = path_vector["hypoPerChr_rds_path"],
    dmrPerChr_plot = path_vector["hypoPerChr_plot_path"],
    dmr_annotated_path = path_vector["hypo_annotated_rds_path"],
    gene_part_annot_plot = path_vector["hypo_gene_part_annotation_plot"],
    tss_association = path_vector["hypo_tss_association_table"])

  diff_meth_analysis(diffMethObj, direction="hyper",
    dmr_path = path_vector["hyper_rds_path"],
    dmrPerChr_path = path_vector["hyperPerChr_rds_path"],
    dmrPerChr_plot = path_vector["hyperPerChr_plot_path"],
    dmr_annotated_path = path_vector["hyper_annotated_rds_path"],
    gene_part_annot_plot = path_vector["hyper_gene_part_annotation_plot"],
    tss_association = path_vector["hyper_tss_association_table"])

}

export_bed <- function(paths) {
    for (p in paths) {
        bnp <- basename(p)
        obj <- readRDS(p)
        gr <- as(obj, "GRanges")
        bed_path <- file.path(bed_folder, sub("rds", "bed", bnp))
        export(gr, con = bed_path, format = "BED")
        message("%s exported to %s", bnp, bed_path)
    }
}




########
# MAIN #
########

##--- load sample sheet, get sample_ids, input paths and treatment vector
sample_sheet <- read.csv(samples_sheet_path, header = TRUE, sep=";", stringsAsFactor=FALSE)
sample.id    <- sample_sheet$sample_name
input_paths  <- file.path(alignments_folder, paste0(sample.id, mate1_pattern, "_val_1_bismark_bt2_pe.CpG_report.txt.gz"))
treatment    <- as.numeric(sample_sheet$treatment)
treatment_levels <- unique(treatment)

if (length(treatment_levels) > 2) {
    message("Detected more than two experimental treatments. Making separate MethylKitDB objects for each comparison.")

    for (i in seq(1,max(unique(treatment_levels)))) {

        cur_treatment <- treatment[treatment %in% c(0, i)]
        cur_treatment <- as.numeric(cur_treatment != 0)
        cur_sample_id <- subset(sample_sheet, treatment %in% c(0,i))$sample_name

        cur_methylRawObj_rds_path    <- add_number_to_path(methylRawObj_rds_path, i)
        cur_methylMergedObj_rds_path <- add_number_to_path(methylMergedObj_rds_path, i)
        cur_correlogram_path         <- add_number_to_path(correlogram_path, i)
        cur_clustering_path          <- add_number_to_path(clustering_path, i)
        cur_diffMeth_rds_path        <- add_number_to_path(diffMeth_rds_path, i)
        cur_tileCounts_rds_path      <- add_number_to_path(tileCounts_rds_path, i)
        cur_diffMethTiles_rds_path   <- add_number_to_path(diffMethTiles_rds_path, i)

        cur_dmc_paths <- sapply(seq_along(dmc_paths), function(ii) add_number_to_path(dmc_paths[ii], i))
        cur_dmr_paths <- sapply(seq_along(dmr_paths), function(ii) add_number_to_path(dmr_paths[ii], i))


        message(sprintf("Using %s as treatment and %s as control.",
          paste(cur_sample_id[cur_treatment == 1], collapse=", "),
          paste(cur_sample_id[cur_treatment == 0], collapse=", ")))

        cur_input_paths <- file.path(alignments_folder, paste0(cur_sample_id, mate1_pattern, "_val_1_bismark_bt2_pe.CpG_report.txt.gz"))

        methylRawObj <- make_methylkitdb_object(cur_input_paths,
          cur_sample_id, cur_treatment,
          methylRawObj_rds_path = cur_methylRawObj_rds_path,
          .methylkitdb_folder = methylkitdb_folder)

        methylMergedObj <- merge_methylRawObj(methylRawObj,
          .methylMergedObj_rds_path = cur_methylMergedObj_rds_path,
          .correlogram_path = cur_correlogram_path,
          .clustering_path = cur_clustering_path)

        diffMeth <- compute_diffMeth(methylMergedObj, diffMethObj_rds_path = cur_diffMeth_rds_path)


        #########
        ## DMC ##
        #########

        do_diff_meth_analysis(diffMeth, cur_dmc_paths)
        export_bed(cur_dmc_paths)

        #########
        ## DMR ##
        #########

        tiles <- compute_tiles(methylMergedObj, .tileCounts_rds_path = cur_tileCounts_rds_path)

        diffMethTiles <- compute_diffMeth(tiles, diffMethObj_rds_path = cur_diffMethTiles_rds_path)

        do_diff_meth_analysis(diffMethTiles, cur_dmr_paths)
        export_bed(cur_dmr_paths)

    }

} else {

    methylRawObj <- make_methylkitdb_object(input_paths, sample.id, treatment, methylRawObj_rds_path)
    methylMergedObj <- merge_methylRawObj(methylRawObj)

    diffMeth <- compute_diffMeth(methylMergedObj)

    #########
    ## DMC ##
    #########

    do_diff_meth_analysis(diffMeth, dmc_paths)
    export_bed(cur_dmc_paths)

    #########
    ## DMR ##
    #########

    tiles <- compute_tiles(methylMergedObj)

    diffMethTiles <- compute_diffMeth(tiles, diffMethObj_rds_path = diffMethTiles_rds_path)

    do_diff_meth_analysis(diffMethTiles, dmr_paths)
    export_bed(cur_dmr_paths)

}

message("all done.")
