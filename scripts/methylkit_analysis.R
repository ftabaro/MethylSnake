sink(stdout(), type = "message")
message("#################################\n## Starting methylKit analysis ##\n#################################")

suppressPackageStartupMessages(library(methylKit))
suppressPackageStartupMessages(library(genomation))
suppressPackageStartupMessages(library(GenomicRanges))



######################
## READ CONFIG FILE ##
######################

# shared folders
methylkitdb_folder <- snakemake@config[["methylkitdb_folder"]]
alignments_folder  <- snakemake@config[["alignments_folder"]]

rdata_folder    <- snakemake@config[["rdata_folder"]]
pictures_folder <- snakemake@config[["pictures_folder"]]
tables_folder   <- snakemake@config[["tables_folder"]]
message(sprintf("Folder summary:\n\tmethylkitdb_folder: %s\n\talignments_folder: %s\n\trdata_folder: %s\n\tpictures_folder: %s\n\ttables_folder: %s",methylkitdb_folder, alignments_folder, rdata_folder, pictures_folder, tables_folder))

# dmr paramters
window_size <- as.numeric(snakemake@config[["dmr_window_size"]])
step_size   <- as.numeric(snakemake@config[["dmr_step_size"]])
diff        <- as.numeric(snakemake@config[["dmr_difference"]])
qvalue      <- as.numeric(snakemake@config[["dmr_qvalue"]])
message(sprintf("DMR calling parameters: wd=%d (%s) - ss=%d  (%s) - diff = %d (%s) - q-value=%.6f (%s)", window_size, class(window_size), step_size, class(step_size), diff, class(diff), qvalue, class(qvalue)))

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

dmr_rds_path      <- file.path(dmr_rdata_folder, "dmr.rds")
hyperDmr_rds_path <- file.path(dmr_rdata_folder, "hyper.rds")
hypoDmr_rds_path  <- file.path(dmr_rdata_folder, "hypo.rds")

dmrPerChr_rds_path      <- file.path(dmr_rdata_folder, "dmrPerChromosome.rds")
hypoDmrPerChr_rds_path  <- file.path(dmr_rdata_folder, "hypoDmrPerChromosome.rds")
hyperDmrPerChr_rds_path <- file.path(dmr_rdata_folder, "hyperDmrPerChromosome.rds")

dmr_annotated_rds_path      <- file.path(dmr_rdata_folder, "dmr_annotated.rds")
hypoDmr_annotated_rds_path  <- file.path(dmr_rdata_folder, "hypo_dmr_annotated.rds")
hyperDmr_annotated_rds_path <- file.path(dmr_rdata_folder, "hyper_dmr_annotated.rds")

dmc_rds_path      <- file.path(dmc_rdata_folder, "dmc.rds")
hyperDmc_rds_path <- file.path(dmc_rdata_folder, "hyper.rds")
hypoDmc_rds_path  <- file.path(dmc_rdata_folder, "hypo.rds")

dmcPerChr_rds_path      <- file.path(dmc_rdata_folder, "dmcPerChromosome.rds")
hypoDmcPerChr_rds_path  <- file.path(dmc_rdata_folder, "hypoDmcPerChromosome.rds")
hyperDmcPerChr_rds_path <- file.path(dmc_rdata_folder, "hyperDmcPerChromosome.rds")

dmc_annotated_rds_path      <- file.path(dmc_rdata_folder, "dmc_annotated.rds")
hypoDmc_annotated_rds_path  <- file.path(dmc_rdata_folder, "hypoDmc_annotated.rds")
hyperDmc_annotated_rds_path <- file.path(dmc_rdata_folder, "hyperDmc_annotated.rds")


# plot paths
correlogram_path              <- file.path(pictures_folder, "samples_correlation.pdf")
clustering_path               <- file.path(pictures_folder, "samples_clustering_pca.pdf")

dmr_gene_part_annotation_plot       <- file.path(dmr_pictures_folder, "annotation.pdf")
hypoDmr_gene_part_annotation_plot  <- file.path(dmr_pictures_folder, "hypo_annotation.pdf")
hyperDmr_gene_part_annotation_plot <- file.path(dmr_pictures_folder, "hyper_annotation.pdf")

dmrPerChr_plot_path       <- file.path(dmr_pictures_folder, "dmrPerChromosome.pdf")
hypoDmrPerChr_plot_path  <- file.path(dmr_pictures_folder, "hypoDmrPerChromosome.pdf")
hyperDmrPerChr_plot_path <- file.path(dmr_pictures_folder, "hyperDmrPerChromosome.pdf")

dmc_gene_part_annotation_plot       <- file.path(dmc_pictures_folder, "annotation.pdf")
hypoDmc_gene_part_annotation_plot  <- file.path(dmc_pictures_folder, "hypo_annotation.pdf")
hyperDmc_gene_part_annotation_plot <- file.path(dmc_pictures_folder, "hyper_annotation.pdf")

dmcPerChr_plot_path       <- file.path(dmc_pictures_folder, "dmcPerChromosome.pdf")
hypoDmcPerChr_plot_path  <- file.path(dmc_pictures_folder, "hypoDmcPerChromosome.pdf")
hyperDmcPerChr_plot_path <- file.path(dmc_pictures_folder, "hyperDmcPerChromosome.pdf")

# table paths
dmr_tss_association_table      <- file.path(dmr_tables_folder, "dmr_tss_association.csv")
hypoDmr_tss_association_table  <- file.path(dmr_tables_folder, "hypoDmr_tss_association.csv")
hyperDmr_tss_association_table <- file.path(dmr_tables_folder, "hyperDmr_tss_association.csv")

dmc_tss_association_table      <- file.path(dmc_tables_folder, "dmc_tss_association.csv")
hypoDmc_tss_association_table  <- file.path(dmc_tables_folder, "hypoDmc_tss_association.csv")
hyperDmc_tss_association_table <- file.path(dmc_tables_folder, "hyperDmc_tss_association.csv")



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
  message("MethylKit database folder is: ", methylkitdb_folder)

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

    meth <- unite(methylRawObj)
    saveRDS(meth, .methylMergedObj_rds_path)

    pdf(.correlogram_path, height=30, width=30)
    getCorrelation(meth, plot = TRUE)
    dev.off()

    pdf(.clustering_path, height=7, width=21)
    par(pty="s", mfrow=c(1,3))
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
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



call_dmr <- function(diffMeth, .dmr_rds_path = dmr_rds_path, .dmrPerChr_rds_path = dmrPerChr_rds_path, .dmrPerChr_plot_path = dmrPerChr_plot_path, direction=NULL){

    if (is.null(direction)){
      direction <- "all"
    }

    dmr <- getMethylDiff(diffMeth,
      difference=diff,
      qvalue=qvalue,
      type=direction)
    saveRDS(dmr, .dmr_rds_path)

    dmrPerChr <- diffMethPerChr(diffMeth, plot = FALSE,
      qvalue.cutoff=qvalue, meth.cutoff=diff )
    saveRDS(dmrPerChr, .dmrPerChr_rds_path)

    chroms <- levels(seqnames(as(diffMeth, "GRanges")))
    exclude <- as.factor(grep("_", chroms, value = TRUE))

    pdf(.dmrPerChr_plot_path)
    diffMethPerChr(diffMeth, plot = TRUE,
      qvalue.cutoff=qvalue,
      meth.cutoff=diff,
      exclude = exclude)
    dev.off()

    message("DMR computation complete")
    return(dmr)
}



annotate <- function(diffMeth, .dmr_annotated_rds_path = dmr_annotated_rds_path, .gene_part_annotation_plot = gene_part_annotation_plot, .tss_association_table = tss_association_table){

    dmr_gr <- as(diffMeth, "GRanges")

    seqlevels(dmr_gr, pruning.mode="coarse") <- seqlevels(gene.obj)

    annotated <- annotateWithGeneParts(dmr_gr, gene.obj)
    saveRDS(annotated, .dmr_annotated_rds_path)

    pdf(.gene_part_annotation_plot)
    plotTargetAnnotation(annotated,
      precedence = TRUE,
      main="Annotation of DMRs")
    dev.off()

    tbl <- getAssociationWithTSS(annotated)
    write.csv(tbl, file=.tss_association_table)

    message("Annotation complete.")
}



add_number_to_path <- function(path, i) {
  sub("(.*)\\.(.*)$", sprintf("\\1\\.%d.\\2", i), path)
}



diff_meth_analysis <- function (methDiffObj,
  dmr_path, dmrPerChr_path, dmrPerChr_plot,
  dmr_annotated_path, gene_part_annot_plot, tss_association,
  direction = "all") {

  dmr <- call_dmr(methDiffObj,
    direction = direction,
    .dmr_rds_path = dmr_path,
    .dmrPerChr_rds_path = dmrPerChr_path,
    .dmrPerChr_plot_path = dmrPerChr_plot)


  annotate(dmr,
    .dmr_annotated_rds_path = dmr_annotated_path,
    .gene_part_annotation_plot = gene_part_annot_plot,
    .tss_association_table = tss_association)
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

        cur_treatment <- treatment[treatment %in% c(0,i)]
        cur_treatment <- as.numeric(cur_treatment != 0)
        cur_sample_id <- subset(sample_sheet, treatment %in% c(0,i))$sample_name

        message(sprintf("Using %s as treatment and %s as control.",
          paste(cur_sample_id[cur_treatment == 1], collapse=", "),
          paste(cur_sample_id[cur_treatment == 0], collapse=", ")))

        cur_input_paths <- file.path(alignments_folder, paste0(cur_sample_id, mate1_pattern, "_val_1_bismark_bt2_pe.CpG_report.txt.gz"))

        methylRawObj <- make_methylkitdb_object(cur_input_paths,
          cur_sample_id,
          cur_treatment,
          methylRawObj_rds_path = add_number_to_path(methylRawObj_rds_path, i),
          .methylkitdb_folder = methylkitdb_folder)

        methylMergedObj <- merge_methylRawObj(methylRawObj,
          .methylMergedObj_rds_path = add_number_to_path(methylMergedObj_rds_path, i),
          .correlogram_path = add_number_to_path(correlogram_path, i),
          .clustering_path = add_number_to_path(clustering_path, i))

        diffMeth <- compute_diffMeth(methylMergedObj,
          diffMethObj_rds_path = add_number_to_path(diffMeth_rds_path, i))

        #########
        ## DMC ##
        #########

        diff_meth_analysis(diffMeth,
          dmr_path = add_number_to_path(dmc_rds_path, i),
          dmrPerChr_path = add_number_to_path(dmcPerChr_rds_path, i),
          dmrPerChr_plot = add_number_to_path(dmcPerChr_plot_path, i),
          dmr_annotated_path = add_number_to_path(dmc_annotated_rds_path, i),
          gene_part_annot_plot = add_number_to_path(dmc_gene_part_annotation_plot, i),
          tss_association = add_number_to_path(dmc_tss_association_table, i))


        ##############
        ## Hypo DMC ##
        ##############

        diff_meth_analysis(diffMeth,
          direction="hypo",
          dmr_path = add_number_to_path(hypoDmc_rds_path, i),
          dmrPerChr_path = add_number_to_path(hypoDmcPerChr_rds_path, i),
          dmrPerChr_plot = add_number_to_path(hypoDmcPerChr_plot_path, i),
          dmr_annotated_path = add_number_to_path(hypoDmc_annotated_rds_path, i),
          gene_part_annot_plot = add_number_to_path(hypoDmc_gene_part_annotation_plot, i),
          tss_association = add_number_to_path(hypoDmc_tss_association_table, i))


        ###############
        ## Hyper DMC ##
        ###############

        diff_meth_analysis(diffMeth,
          direction="hyper",
          dmr_path = add_number_to_path(hyperDmc_rds_path, i),
          dmrPerChr_path = add_number_to_path(hyperDmcPerChr_rds_path, i),
          dmrPerChr_plot = add_number_to_path(hyperDmcPerChr_plot_path, i),
          dmr_annotated_path = add_number_to_path(hyperDmc_annotated_rds_path, i),
          gene_part_annot_plot = add_number_to_path(hyperDmc_gene_part_annotation_plot, i),
          tss_association = add_number_to_path(hyperDmc_tss_association_table, i))



        #########
        ## DMR ##
        #########

        tiles <- compute_tiles(methylMergedObj,
          .tileCounts_rds_path = add_number_to_path(tileCounts_rds_path, i))

        diffMethTiles <- compute_diffMeth(tiles,
          diffMethObj_rds_path = add_number_to_path(diffMethTiles_rds_path, i))

        diff_meth_analysis(diffMethTiles,
          dmr_path = add_number_to_path(dmr_rds_path, i),
          dmrPerChr_path = add_number_to_path(dmrPerChr_rds_path, i),
          dmrPerChr_plot = add_number_to_path(dmrPerChr_plot_path, i),
          dmr_annotated_path = add_number_to_path(dmr_annotated_rds_path, i),
          gene_part_annot_plot = add_number_to_path(dmr_gene_part_annotation_plot, i),
          tss_association = add_number_to_path(dmr_tss_association_table, i))

        ##############
        ## Hypo DMR ##
        ##############

        diff_meth_analysis(diffMethTiles,
          direction="hypo",
          dmr_path = add_number_to_path(hypoDmr_rds_path, i),
          dmrPerChr_path = add_number_to_path(hypoDmrPerChr_rds_path, i),
          dmrPerChr_plot = add_number_to_path(hypoDmrPerChr_plot_path, i),
          dmr_annotated_path = add_number_to_path(hypoDmr_annotated_rds_path, i),
          gene_part_annot_plot = add_number_to_path(hypoDmr_gene_part_annotation_plot, i),
          tss_association = add_number_to_path(hypoDmr_tss_association_table, i))

        ###############
        ## Hyper DMR ##
        ###############

        diff_meth_analysis(diffMethTiles,
          direction="hyper",
          dmr_path = add_number_to_path(hyperDmr_rds_path, i),
          dmrPerChr_path = add_number_to_path(hyperDmrPerChr_rds_path, i),
          dmrPerChr_plot = add_number_to_path(hyperDmrPerChr_plot_path, i),
          dmr_annotated_path = add_number_to_path(hyperDmr_annotated_rds_path, i),
          gene_part_annot_plot = add_number_to_path(hyperDmr_gene_part_annotation_plot, i),
          tss_association = add_number_to_path(hyperDmr_tss_association_table, i))

    }

} else {

    methylRawObj <- make_methylkitdb_object(input_paths, sample.id, treatment, methylRawObj_rds_path)
    methylMergedObj <- merge_methylRawObj(methylRawObj)

    diffMeth <- compute_diffMeth(methylMergedObj,
      diffMethObj_rds_path = add_number_to_path(diffMeth_rds_path, i))

    #########
    ## DMC ##
    #########

    diff_meth_analysis(diffMeth,
      dmr_path = dmc_rds_path,
      dmrPerChr_path = dmcPerChr_rds_path,
      dmrPerChr_plot = dmcPerChr_plot_path,
      dmr_annotated_path = dmc_annotated_rds_path,
      gene_part_annot_plot = dmc_gene_part_annotation_plot,
      tss_association = dmc_tss_association_table)

    ##############
    ## Hypo DMC ##
    ##############

    diff_meth_analysis(diffMeth,
      direction="hypo",
      dmr_path = hypoDmc_rds_path,
      dmrPerChr_path = hypoDmcPerChr_rds_path,
      dmrPerChr_plot = hypoDmcPerChr_plot_path,
      dmr_annotated_path = hypoDmc_annotated_rds_path,
      gene_part_annot_plot = hypoDmc_gene_part_annotation_plot,
      tss_association = hypoDmc_tss_association_table)

    ###############
    ## Hyper DMC ##
    ###############

    diff_meth_analysis(diffMeth,
      direction="hyper",
      dmr_path = hyperDmc_rds_path,
      dmrPerChr_path = hyperDmcPerChr_rds_path,
      dmrPerChr_plot = hyperDmcPerChr_plot_path,
      dmr_annotated_path = hyperDmc_annotated_rds_path,
      gene_part_annot_plot = hyperDmc_gene_part_annotation_plot,
      tss_association = hyperDmc_tss_association_table)



    #########
    ## DMR ##
    #########

    tiles <- compute_tiles(methylMergedObj)

    diffMethTiles <- compute_diffMeth(tiles,
      diffMethObj_rds_path = add_number_to_path(diffMethTiles_rds_path, i))

    diff_meth_analysis(tiles,
      dmr_path = dmr_rds_path,
      dmrPerChr_path = dmrPerChr_rds_path,
      dmrPerChr_plot = dmrPerChr_plot_path,
      dmr_annotated_path = dmr_annotated_rds_path,
      gene_part_annot_plot = dmr_gene_part_annotation_plot,
      tss_association = dmr_tss_association_table)

    ##############
    ## Hypo DMR ##
    ##############

    diff_meth_analysis(tiles,
      direction="hypo",
      dmr_path = hypoDmr_rds_path,
      dmrPerChr_path = hypoDmrPerChr_rds_path,
      dmrPerChr_plot = hypoDmrPerChr_plot_path,
      dmr_annotated_path = hypoDmr_annotated_rds_path,
      gene_part_annot_plot = hypoDmr_gene_part_annotation_plot,
      tss_association = hypoDmr_tss_association_table)

    ###############
    ## Hyper DMR ##
    ###############

    diff_meth_analysis(tiles,
      direction="hyper",
      dmr_path = hyperDmr_rds_path,
      dmrPerChr_path = hyperDmrPerChr_rds_path,
      dmrPerChr_plot = hyperDmrPerChr_plot_path,
      dmr_annotated_path = hyperDmr_annotated_rds_path,
      gene_part_annot_plot = hyperDmr_gene_part_annotation_plot,
      tss_association = hyperDmr_tss_association_table)

}

message("all done.")
