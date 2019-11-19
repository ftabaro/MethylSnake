import os

singularity: "docker://continuumio/miniconda3:4.7.10"
#localrules: bismark2report, bismark2summary, bam2nuc, filter_incomplete_conversions
localrules: bismark2report, bismark2summary

mate1_root = "{sample}" + config["mate1_pattern"]
mate2_root = "{sample}" + config["mate2_pattern"]

rule all:
    message: "all done!"
    input:
      expand(os.path.join(config["methylkitdb_folder"], "{sample}.txt.bgz.tbi"), sample=config["samples"]),
      expand(os.path.join(config["reports_folder"], mate1_root + "_val_1_bismark_bt2_PE_report.html"), sample=config["samples"]),
      os.path.join(config["alignments_folder"], "bismark_summary_report.html"),
      expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nonCG_filtered.bam"), sample=config["samples"]),
      expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nonCG_removed_seqs.bam"), sample=config["samples"]),
      expand(os.path.join(config["pictures_folder"], "{sample}_preliminary_stats.pdf"), sample=config["samples"]),
      os.path.join(config["rdata_folder"], "methylMergedObj.rds"),
      os.path.join(config["rdata_folder"], "dmr_annotated.rds")


# rule fastqc:
# rule multiqc:
# rule bismark_genome_prep:

rule trim:
    message: "Performing reads trimming..."
    input:
      [ os.path.join(config["reads_folder"], mate1_root + config["fastq_extension"]), os.path.join(config["reads_folder"], mate2_root + config["fastq_extension"]) ]
    output:
      os.path.join(config["trimmed_folder"], mate1_root  + "_val_1.fq.gz"),
      os.path.join(config["trimmed_folder"], mate1_root  + ".fastq.gz_trimming_report.txt"),
      os.path.join(config["trimmed_folder"], mate2_root  + "_val_2.fq.gz"),
      os.path.join(config["trimmed_folder"], mate2_root  + ".fastq.gz_trimming_report.txt"),
    conda:
      os.path.join(config["environments_folder"], "rrbs.yaml")
    params:
      quality_filter_value="22",
    log:
      os.path.join(config["log_folder"], "trim_galore", "{sample}.log")
    benchmark:
      os.path.join(config["log_folder"], "trim_galore", "{sample}.benchmark.log")
    threads: 8
    shell:
      """
      OUTPUT=$(dirname {output[0]})

      trim_galore --quality {params.quality_filter_value} --phred33 --output_dir $OUTPUT --gzip --rrbs --fastqc --fastqc_args "-o \"$OUTPUT\"" --paired --cores {threads} {input}
      """


rule bismark_align:
    message: "Performing alignment..."
    input:
      os.path.join(config["trimmed_folder"], mate1_root + "_val_1.fq.gz"),
      os.path.join(config["trimmed_folder"], mate2_root + "_val_2.fq.gz")
    output:
      os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.bam"),
      os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_PE_report.txt"),
      os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nucleotide_stats.txt")
    params:
      config["tmp_folder"]
    log:
      os.path.join(config["log_folder"], "bismark", "{sample}.log")
    conda: os.path.join(config["environments_folder"], "rrbs.yaml")
    threads: 4
    benchmark:
      os.path.join(config["log_folder"], "bismark", "{sample}.benchmark.log")
    shell:
      """
      bismark --phred33-quals --bowtie2 -p {threads} --genome {config[bismark_index_path]} --unmapped --ambiguous --ambig_bam --nucleotide_coverage --output_dir {config[alignments_folder]} --fastq --temp_dir {params[0]} -1 {input[0]} -2 {input[1]}
      """

# rule samtools_sort:
#     input:
#       os.path.join(config["alignments_folder"], "{sample}_1_val_1_bismark_bt2_pe.bam")
#     output:
#       os.path.join(config["alignments_folder"], "{sample}_1_val_1_bismark_bt2_pe_sorted.bam")
#     params:
#       tmp_folder=config["tmp_folder"]
#     threads: 8
#     conda: os.path.join(config["environments_folder"], "rrbs.yaml")
#     shell:
#       """
#       samtools sort -o {output} -T {params.tmp_folder} -@ {threads} {input}
#       """

# rule samtools_index:
#   input:
#     os.path.join(config["alignments_folder"], "{sample}_1_val_1_bismark_bt2_pe_sorted.bam")
#   output:
#     os.path.join(config["alignments_folder"], "{sample}_1_val_1_bismark_bt2_pe_sorted.bam.bai")
#   conda: os.path.join(config["environments_folder"], "rrbs.yaml")
#   shell:
#     """
#     samtools index {input}
#     """


rule methylation_extractor:
    message: "Performing methylation extraction..."
    input:
      os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.bam")
    output:
      os.path.join(config["alignments_folder"], "CpG_context_"+ mate1_root +"_val_1_bismark_bt2_pe.txt.gz"),
      os.path.join(config["alignments_folder"], "CHG_context_"+ mate1_root +"_val_1_bismark_bt2_pe.txt.gz"),
      os.path.join(config["alignments_folder"], "CHH_context_"+ mate1_root +"_val_1_bismark_bt2_pe.txt.gz"),
      os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.CpG_report.txt.gz")
    log:
      os.path.join(config["log_folder"], "bismark_methylation_extractor", "{sample}.log")
    conda:
      os.path.join(config["environments_folder"], "rrbs.yaml")
    benchmark:
      os.path.join(config["log_folder"], "bismark_methylation_extractor", "{sample}.benchmark.log")
    threads: 9
    shell:
      """
      bismark_methylation_extractor --paired-end --comprehensive --gzip --output {config[alignments_folder]} --multicore $(( {threads} / 3 )) --bedGraph --remove_spaces --buffer_size 80% --cytosine_report --genome_folder $(dirname {config[genome_fasta]}) --ignore_r2 2 {input[0]}
      """


rule bismark2report:
  message: "Generating reports... "
  input:
      expand(os.path.join(config["alignments_folder"], "CpG_context_"+ mate1_root + "_val_1_bismark_bt2_pe.txt.gz"), sample=config["samples"]),
      expand(os.path.join(config["alignments_folder"], "CHG_context_"+ mate1_root + "_val_1_bismark_bt2_pe.txt.gz"), sample=config["samples"]),
      expand(os.path.join(config["alignments_folder"], "CHH_context_"+ mate1_root + "_val_1_bismark_bt2_pe.txt.gz"), sample=config["samples"]),
      expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nucleotide_stats.txt"), sample=config["samples"])
  output:
      os.path.join(config["reports_folder"], mate1_root + "_val_1_bismark_bt2_PE_report.html")
  log:
      os.path.join(config["log_folder"], "bismark2report", "{sample}.log")
  conda:
      os.path.join(config["environments_folder"], "rrbs.yaml")
  benchmark:
      os.path.join(config["log_folder"], "bismark2report", "{sample}.benchmark.log")
  threads: 1
  shell:
      """
      cd {config[alignments_folder]}
      bismark2report --dir {config[reports_folder]}
      """


rule bismark2summary:
  message: "Generating summary... "
  input:
    expand(os.path.join(config["alignments_folder"], "CpG_context_"+ mate1_root + "_val_1_bismark_bt2_pe.txt.gz"), sample=config["samples"]),
    expand(os.path.join(config["alignments_folder"], "CHG_context_"+ mate1_root + "_val_1_bismark_bt2_pe.txt.gz"), sample=config["samples"]),
    expand(os.path.join(config["alignments_folder"], "CHH_context_"+ mate1_root + "_val_1_bismark_bt2_pe.txt.gz"), sample=config["samples"]),
  output:
    os.path.join(config["alignments_folder"], "bismark_summary_report.html")
  log:
    os.path.join(config["log_folder"], "bismark2summary", "bismark_summary_report.log")
  conda:
    os.path.join(config["environments_folder"], "rrbs.yaml")
  benchmark:
      os.path.join(config["log_folder"], "bismark2summary", "bismark_summary_report.benchmark.log")
  threads: 1
  shell:
      """
      cd {config[alignments_folder]}
      bismark2summary
      """


rule filter_incomplete_conversions:
  message: "Filtering incomplete conversions..."
  input:
    os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.bam")
  output:
    os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nonCG_filtered.bam"),
    os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nonCG_removed_seqs.bam")
  log:
    os.path.join(config["log_folder"], "filter_non_conversion", "{sample}.log")
  conda:
    os.path.join(config["environments_folder"], "rrbs.yaml")
  benchmark:
    os.path.join(config["log_folder"], "filter_non_conversion", "{sample}.benchmark.log")
  threads: 1
  shell:
    """
    filter_non_conversion --paired {input[0]}
    """


#rule download_ensembl:
#  message:
#    "Downloading Ensembl annotations..."
#  output:
#    os.path.join(config["ensembl_root"], config["ensembl_version"], "gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz")
#  params:
#    ensembl_version=config["ensembl_version"]
#  conda:
#    os.path.join(config["environments_folder"], "rrbs.yaml")
#  log:
#    os.path.join(config["log_folder"], "download_ensembl.log")
#  threads: 1
#  shell:
#    """
#    rsync -avz --relative -L rsync://ftp.ensembl.org/pub/{params.ensembl_version} {output}
#    """


rule convert_gtf_to_bed12:
  message: "Converting Ensembl GTF to BED"
  input:
#    os.path.join(config["ensembl_root"], config["ensembl_version"], "gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz")
    config["annotation_file"]
  output:
    re.sub("gtf", "bed12", config["annotation_file"])
#    os.path.join(config["ensembl_root"], config["ensembl_version"], "gtf/homo_sapiens/Homo_sapiens.GRCh37.87.bed12.gz")
  conda:
    os.path.join(config["environments_folder"], "rrbs.yaml")
  params:
    tmp_dir=config["tmp_folder"]
  shell:
    """
    cd {params.tmp_dir}
    gunzip -c {input} > temp.gtf;
    gtfToGenePred temp.gtf temp.genePred;
    genePredToBed temp.genePred temp.bed;
    gzip -c temp.bed > {output};
    rm temp.gtf temp.genePred temp.bed
    """


rule make_methylkit_db:
  message:
    "Generating MethylKit database..."
  input:
    expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.CpG_report.txt.gz"), sample=config["samples"]),
    sample_sheet=config["sample_sheet"]
  output:
    methylRawObj=os.path.join(config["rdata_folder"], "methylRawObj.rds"),
    bgzFiles=expand(os.path.join(config["methylkitdb_folder"], "{sample}.txt.bgz.tbi"), sample=config["samples"]),
    plots=expand(os.path.join(config["pictures_folder"], "{sample}_preliminary_stats.pdf"), sample=config["samples"])
  params:
    methylkitdb_folder=config["methylkitdb_folder"],
    pictures_folder=config["pictures_folder"]
  conda:
    os.path.join(config["environments_folder"], "methylkit.yaml")
  script:
    "scripts/methylkit/make_db.R"


rule methylkit_merge_samples:
  message: "Merging samples..."
  input:
    methylRawObj=os.path.join(config["rdata_folder"], "methylRawObj.rds"),
    dbFiles=expand(os.path.join(config["methylkitdb_folder"], "{sample}.txt.bgz.tbi"), sample=config["samples"])
  output:
    methylMergedObj=os.path.join(config["rdata_folder"], "methylMergedObj.rds"),
    correlogram=os.path.join(config["pictures_folder"], "samples_correlation.pdf"),
    clustering=os.path.join(config["pictures_folder"], "samples_clustering_pca.pdf")
  conda:
    os.path.join(config["environments_folder"], "methylkit.yaml")
  script:
    "scripts/methylkit/merge_samples.R"


rule methylkit_detect_dmr:
  message: "Performing tiling windows analysis to detect DMR..."
  input:
    methylMergedObj=os.path.join(config["rdata_folder"], "methylMergedObj.rds")
  output:
    tileCounts=os.path.join(config["rdata_folder"], "tilesCounts.rds"),
    dmr=os.path.join(config["rdata_folder"], "dmr.rds"),
    dmrPerChromosome=os.path.join(config["rdata_folder"], "dmrPerChromosome.rds"),
    dmrHyper=os.path.join(config["rdata_folder"], "dmrHyper.rds"),
    dmrHypo=os.path.join(config["rdata_folder"], "dmrHypo.rds"),
    dmrPerChromosome_plot=os.path.join(config["pictures_folder"], "dmrPerChromosome.pdf")
  threads: 6
  params:
    window_size=config["dmr_window_size"],
    step_size=config["dmr_step_size"],
    difference=config["dmr_difference"],
    qvalue=config["dmr_qvalue"]
  conda:
    os.path.join(config["environments_folder"], "methylkit.yaml")
  script:
    "scripts/methylkit/detect_dmr.R"

rule annotate_dmr:
  input:
    dmr=os.path.join(config["rdata_folder"], "dmr.rds"),
#    annotation=os.path.join(config["ensembl_root"], config["ensembl_version"], "gtf/homo_sapiens/Homo_sapiens.GRCh37.87.bed12.gz")
    annotation=re.sub("gtf", "bed12", config["annotation_file"])
  output:
    dmr_annotated=os.path.join(config["rdata_folder"], "dmr_annotated.rds"),
    gene_part_annotation_plot=os.path.join(config["pictures_folder"], "dmr_annotation.pdf"),
    tss_association=os.path.join(config["tables_folder"], "tss_association.csv"),
    genomic_features_percentages=os.path.join(config["tables_folder"], "genomic_features_percentages.txt")
  conda:
    os.path.join(config["environments_folder"], "methylkit.yaml")
  script:
    "scripts/methylkit/annotate_dmr.R"
