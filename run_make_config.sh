#!/bin/sh

REFERENCE="/bmt-data/genomics/reference"

python scripts/make_config.py \
    --wd .. \
    --config-path ../config.yaml \
    --genome-path $REFERENCE/hg38.fa \
    --bismark-index-path $REFERENCE/bismark/hg38 \
    --annotation-file $REFERENCE/Ensembl/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz \
    --sample-sheet /lustre/bmt-data/genomics/projects/glioma_cellline_cfrrbs/sample_sheet.csv \
    --tmp-folder /scratch/ft413468/rrbs-sharma \
    --mate1-pattern _1 \
    --mate2-pattern _2 \
    --fastq-extension .fq.gz \
    --log-folder log \
    --reads-folder /lustre/bmt-data/genomics/projects/glioma_cellline_cfrrbs/pooled

