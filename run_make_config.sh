#!/bin/sh

REFERENCE="/bmt-data/genomics/reference"

python scripts/make_config.py --wd .. --config-path ../config.yaml \
  --genome-path $REFERENCE/hg38.fa \
  --bismark-index-path $REFERENCE/bismark/hg38 \
  --annotation-file $REFERENCE/Ensembl/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz \
  --sample-sheet /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/sample_sheet.csv \
  --tmp-folder /scratch/ft413468/rrbs-sharma \
  --mate1-pattern _R1_001 --mate2-pattern _R2_001 --fastq-extension .fastq.gz
