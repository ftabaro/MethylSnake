import os

singularity: config["singularity_container"]

localrules: bismark2report, bismark2summary

mate1_root = "{sample}" + config["mate1_pattern"]
mate2_root = "{sample}" + config["mate2_pattern"]

rule all:
    message: "all done!"
    input:
      expand(os.path.join(config["reports_folder"], mate1_root + "_val_1_bismark_bt2_PE_report.html"), sample=config["samples"]),
      os.path.join(config["alignments_folder"], "bismark_summary_report.html"),
      expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nonCG_filtered.bam"), sample=config["samples"]),
      expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.nonCG_removed_seqs.bam"), sample=config["samples"]),
      os.path.join(config["wd"], "methylkit_analysis.done"),
      expand(os.path.join(config["fastqc_folder"], mate1_root + "_fastqc.zip"), sample=config["samples"]),
      expand(os.path.join(config["fastqc_folder"], mate1_root + "_fastqc.html"), sample=config["samples"]),
      expand(os.path.join(config["fastqc_folder"], mate2_root + "_fastqc.zip"), sample=config["samples"]),
      expand(os.path.join(config["fastqc_folder"], mate2_root + "_fastqc.html"), sample=config["samples"]),

 
rule fastqc:
    message: "Running FastQC..."
    input:
        os.path.join(config["reads_folder"], mate1_root + config["fastq_extension"]),
        os.path.join(config["reads_folder"], mate2_root + config["fastq_extension"]),
    output:
        os.path.join(config["fastqc_folder"], mate1_root + "_fastqc.zip"),
        os.path.join(config["fastqc_folder"], mate1_root + "_fastqc.html"),
        os.path.join(config["fastqc_folder"], mate2_root + "_fastqc.zip"),
        os.path.join(config["fastqc_folder"], mate2_root + "_fastqc.html")
    params:
        dir=config["fastqc_folder"]
    threads: 6
    shell:
        """
        fastqc -o {params.dir} -t {threads} --noextract {input[0]}
        fastqc -o {params.dir} -t {threads} --noextract {input[1]}
        """

# rule multiqc:
# rule bismark_genome_prep:

rule trim:
    message: "Performing reads trimming..."
    input:
      [ os.path.join(config["reads_folder"], mate1_root + config["fastq_extension"]), os.path.join(config["reads_folder"], mate2_root + config["fastq_extension"]) ]
    output:
      os.path.join(config["trimmed_folder"], mate1_root  + "_val_1.fq.gz"),
      os.path.join(config["trimmed_folder"], mate2_root  + "_val_2.fq.gz"),
      os.path.join(config["trimmed_folder"], mate1_root  + config["fastq_extension"] + "_trimming_report.txt"),
      os.path.join(config["trimmed_folder"], mate2_root  + config["fastq_extension"] + "_trimming_report.txt"),
    #conda:
    #  os.path.join(config["environments_folder"], "rrbs.yaml")
    params:
      quality_filter_value="22",
    log:
      os.path.join(config["log_folder"], "trim_galore", "{sample}.log")
    benchmark:
      os.path.join(config["log_folder"], "trim_galore", "{sample}.benchmark.log")
    threads: 4
    shell:
      """
      OUTPUT=$(dirname {output[0]})

      trim_galore --quality {params.quality_filter_value} \
      --phred33 --output_dir $OUTPUT \
      --gzip --rrbs --fastqc --fastqc_args "-o \"$OUTPUT\"" \
      --paired --cores {threads} {input}

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
    #conda: os.path.join(config["environments_folder"], "rrbs.yaml")
    threads: 4
    benchmark:
      os.path.join(config["log_folder"], "bismark", "{sample}.benchmark.log")
    shell:
      """
      bismark --phred33-quals --bowtie2 -p {threads} \
        --genome {config[bismark_index_path]} --unmapped --ambiguous \
        --ambig_bam --nucleotide_coverage --output_dir {config[alignments_folder]} \
        --fastq --temp_dir {params[0]} -1 {input[0]} -2 {input[1]}
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
#    conda:
#      os.path.join(config["environments_folder"], "rrbs.yaml")
    benchmark:
      os.path.join(config["log_folder"], "bismark_methylation_extractor", "{sample}.benchmark.log")
    threads: 9
    shell:
      """
      bismark_methylation_extractor --paired-end --comprehensive --gzip \
      --output {config[alignments_folder]} --multicore $(( {threads} / 3 )) \
      --bedGraph --remove_spaces --buffer_size 80% --cytosine_report \
      --genome_folder $(dirname {config[genome_path]}) --ignore_r2 2 {input[0]}
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
    config["annotation_file"]
  output:
    re.sub("gtf.gz$", "bed12", config["annotation_file"])
  params:
    tmp_dir=config["tmp_folder"]
  shell:
    """
    cd {params.tmp_dir}
    gunzip -c {input} > temp.gtf;
    gtfToGenePred temp.gtf temp.genePred;
    genePredToBed temp.genePred {output};
    rm temp.gtf temp.genePred
    """

rule run_methylkit_analysis:
  message:
    "Performing methylKit analysis..."
  input:
    expand(os.path.join(config["alignments_folder"], mate1_root + "_val_1_bismark_bt2_pe.CpG_report.txt.gz"), sample=config["samples"]),
    annotation_file=re.sub("gtf.gz", "bed12", config["annotation_file"]),
    sample_sheet=config["sample_sheet"]
  output:
    touch(os.path.join(config["wd"], "methylkit_analysis.done"))
  script:
    "scripts/methylkit_analysis.R"
