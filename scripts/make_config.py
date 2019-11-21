"""
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

    parser.add_argument("--config-path", required=True)
    parser.add_argument("--wd", required=True)
    parser.add_argument("--genome-path", required=True)
    parser.add_argument("--bismark-index-path", required=True)
    parser.add_argument("--sample-sheet", required=True)

    parser.add_argument("--annotation-file", required=False)
#    parser.add_argument("--ensembl-root", required=False)
#    parser.add_argument("--ensembl-version", required=False)

    parser.add_argument("--environments-folder", default="../environments")
    parser.add_argument("--singularity-container", default="/home/ft413468/containers/methylkit.sif")

    parser.add_argument("--tmp-folder", default="tmp")
    parser.add_argument("--log-folder", default="~/var/log")
    parser.add_argument("--reads-folder", default="reads")
    parser.add_argument("--trimmed-folder", default="trimmed")
    parser.add_argument("--alignments-folder", default="alignments")
    parser.add_argument("--reports-folder", default="reports")
    parser.add_argument("--nucleotide-stats-folder", default="nucleotide_stats")
    parser.add_argument("--methylkitdb-folder", default="methylDB")
    parser.add_argument("--rdata-folder", default="RData")
    parser.add_argument("--pictures-folder", default="pictures")
    parser.add_argument("--tables-folder", default="tables")
    parser.add_argument("--dmr-window-size", default=500)
    parser.add_argument("--dmr-step-size", default=500)
    parser.add_argument("--dmr-difference", default=25)
    parser.add_argument("--dmr-qvalue", default=0.01)

    parser.add_argument("--mate1-pattern", default="_1")
    parser.add_argument("--mate2-pattern", default="_2")
    parser.add_argument("--fastq-extension", default="_2")

    args = parser.parse_args()
    return vars(args)


def validate_args(args):

    for k in ["config_path", "wd", "genome_path", "bismark_index_path", "sample_sheet", "environments_folder"]:
        if not os.path.isabs(args[k]) and args[k] != "-":
            args[k] = os.path.abspath(args[k])

    for k in ["tmp_folder", "log_folder", "reads_folder", "trimmed_folder", "alignments_folder", "reports_folder", "nucleotide_stats_folder", "methylkitdb_folder", "rdata_folder", "pictures_folder", "tables_folder"]:
        if not os.path.isabs(args[k]):
            args[k] = os.path.join(args["wd"], args[k])

    sample_names = []
    with open(args["sample_sheet"], "r") as fin:
        lines = fin.readlines()
        for i in range(1, len(lines)):
            line = lines[i]
            line = line.strip()
            line = re.split(r"\W", line)
            line = [x.replace("\"", "") for x in line]
            line = [x for x in line if len(x) > 0]
            sample_names.append(line[0])
    args["samples"] = sample_names
    return args


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
        args = validate_args(args)
        write_config(args)
