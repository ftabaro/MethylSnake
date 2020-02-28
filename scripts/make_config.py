"""
Author: Francesco Tabaro <francesco.tabaro@gmail.com>

This script will parse its command line arguments and generate a YAML file suitable the RRBS Snakemake pipeline.

It will parse the given sample sheet to generate an array of sample names to use in the alignment phase and

make_config.py \
  --wd /path/to/wd \
  --tmp_folder /path/to/tmp \
  --log_folder /path/to/log \
  --reads-folder /path/to/fastq \
  --trimmed-folder /path/to/trimmed \
  --alignments-folder /path/to/bam \
  --reports-folder /path/to/reports \
  --nucleotide-stats-folder /path/to/nt/stats \
  --genome-path /path/to/genome/fasta \
  --index-folder /path/to/bismark/index \
  --ensembl-root /path/to/ensembl \
  --ensembl-version path/relative/to/ensemlb/root \
  --methylkitdb-folder /path/to/methylkitdb \
  --rdata-folder /parh/to/RData \
  --pictures-folder /path/to/pictures \
  --tables-folder /path/to/tables \
  --tile-size 500 \
  --step-size 500 \
  --dmr-diff 25 \
  --dmr-qvalue 0.01 \
  --sample-sheet /path/to/sample_sheet.csv \
"""

import argparse
import os
import pyaml
import sys
import re


def make_parser():
    parser = argparse.ArgumentParser()

    ## Mandatory options
    parser.add_argument("--config-path",
        required=True,
        help="Desired path for the config file")

    parser.add_argument("--wd",
        required=True,
        help="Working directory")

    parser.add_argument("--genome-path",
        required=True,
        help="Path to a fasta file with genome sequence")

    parser.add_argument("--bismark-index-path",
        required=True,
        help="Path to Bismark index files (only index name required: /path/to/index/folder/hg38*)")

    parser.add_argument("--sample-sheet",
        required=True,
        help="Path to csv file holding samples information")

    parser.add_argument("--annotation-file",
        required=True,
        help="Path to a gzipped GTF file.")

    ## Numeric parameters for DMR detection
    parser.add_argument("--dmr-window-size",
        default=500,
        help="Window size for tiled differential methylation analysis")

    parser.add_argument("--dmr-step-size",
        default=500,
        help="Step size for tiled differential methylation analysis")

    parser.add_argument("--dmr-difference",
        default=25,
        help="Difference in reads coverage threshold for differential methylation analysis")

    parser.add_argument("--dmr-qvalue",
        default=0.01,
        help="Q-value threshold for differential methylation analysis")

    parser.add_argument("--min-per-group",
        default=1,
        help="An integer denoting minimum number of samples per replicate needed to cover a region/base")

    ## Patterns and string parameters
    parser.add_argument("--mate1-pattern",
        default="_1",
        help="Pattern to identify mate 1 in paired sequencing files")

    parser.add_argument("--mate2-pattern",
        default="_2",
        help="Pattern to identify mate 2 in paired sequencing files")

    parser.add_argument("--fastq-extension",
        default=".fq.gz",
        help="File extension of reads fastq files")

    parser.add_argument("--genome-version",
        default="hg38",
        help="Label for genome version used in the analysis")

    ## Optional paths
    parser.add_argument("--singularity-container",
        default="/home/ft413468/containers/methylkit.sif",
        help="Path to a singularity container with all necessary tools installed.")

    parser.add_argument("--tmp-folder",
        default="tmp",
        help="Path to a temporary folder (possibly fast storage)")

    parser.add_argument("--log-folder",
        default="~/var/log",
        help="Path to a log folder")

    parser.add_argument("--reads-folder",
        default="reads",
        help="Path to a folder with reads to be processed")

    parser.add_argument("--trimmed-folder",
        default="trimmed",
        help="Path to folder to write trimmed reads in")

    parser.add_argument("--alignments-folder",
        default="alignments",
        help="Path to a folder to write alignments in")

    parser.add_argument("--reports-folder",
        default="reports",
        help="Path to a folder to write Bismark reports in")

    parser.add_argument("--nucleotide-stats-folder",
        default="nucleotide_stats",
        help="Path to write Bismark nucleotide report files in")

    parser.add_argument("--methylkitdb-folder",
        default="NULL",
        help="Path to write methylKit tabix files in")

    parser.add_argument("--rdata-folder",
        default="RData",
        help="Path to write RDS objects used in methylKit analysis")

    parser.add_argument("--bed-folder",
        default="bed",
        help="Path to write BED files for DMR/DMC coordinates to")

    parser.add_argument("--pictures-folder",
        default="pictures",
        help="Path to a folder to write plots in")

    parser.add_argument("--tables-folder",
        default="tables",
        help="Path to a folder to write tables in")

    parser.add_argument("--fastqc-folder",
        default="fastqc",
        help="FastQC results folder")


    args = parser.parse_args()
    return vars(args)


def test_relative_path(p):
    return not os.path.isabs(p) and not "~" in p


def make_env(args):

    for k in ["config_path", "wd", "genome_path", "bismark_index_path", "sample_sheet"]:
        if test_relative_path(args[k]) and args[k] != "-":
            args[k] = os.path.abspath(args[k])

    for k in ["tmp_folder", "log_folder", "reads_folder", "trimmed_folder", "alignments_folder", "reports_folder", "nucleotide_stats_folder", "methylkitdb_folder", "rdata_folder", "pictures_folder", "tables_folder", "bed_folder", "fastqc_folder"]:

        if args[k] == "NULL":
            continue

        if test_relative_path(args[k]):
            new = os.path.join(args["wd"], args[k])
            print("Setting {}: {} -> {}".format( k, args[k], new))
            args[k] = new

        if not os.path.exists(args[k]):
            os.makedirs(args[k])
            print("Created {}".format(args[k]))
    return args


def parse_samples(sample_sheet_path):
    sample_names = []
    with open(sample_sheet_path, "r") as fin:
        lines = fin.readlines()
        for i in range(1, len(lines)):
            line = lines[i]
            line = line.strip()
            line = re.split(r"\W", line)
            line = [x.replace("\"", "") for x in line]
            line = [x for x in line if len(x) > 0]
            sample_names.append(line[0])
    return sample_names


def write_config(args):
    config_path = args["config_path"]
    del args["config_path"]
    print(config_path)
    if config_path == "-":
        pyaml.dump(args, sys.stdout)
    else:
        fout = open(config_path, "w")
        pyaml.dump(args, fout)
        fout.close()
    return


if __name__ == "__main__":
        args = make_parser()
        args = make_env(args)
        args["samples"] = parse_samples(args["sample_sheet"])
        write_config(args)
