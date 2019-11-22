#!/bin/bash

#if [ ! -d logs ]; then
 # mkdir -p logs/{trim_galore,bismark,bismark_methylation_extractor,bismark2report,bismark2summary,bam2nuc,filter_non_conversion}
#fi

#singularity exec -B /lustre/bmt-data/genomics/projects/glioma_cellline_rrbs ~/containers/snakemake_latest.sif snakemake \
if [ $# -eq 0 ]; then     
    echo "Config file is required."     
    exit 1 
fi

LOG_FOLDER=$(grep log_folder "$1" | awk -F":" '{print $2}')
mkdir -p "$LOG_FOLDER"
mkdir -p $LOG_FOLDER

WD=$(grep wd "$1" | awk -F":" '{print $2}')
GENOME=$(grep genome_path "$1" | awk -F":" '{print $2}' |xargs basename)
BISMARK_IDX=$(grep bismark_index_path "$1" | awk -F":" '{print $2}')
TMP_FOLDER=$(grep tmp_folder "$1" | awk -F":" '{print $2}')

SBATCH="sbatch -J {cluster.jobname} --partition={cluster.partition} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={threads} -e ${LOG_FOLDER%/}/{cluster.log} -o ${LOG_FOLDER%/}/{cluster.log}"

snakemake -r -p -j 100 --configfile "$1" \
    --use-conda \
    --cluster-config ./cluster-config/cluster.json \
    --cluster "$SBATCH" \
    --latency-wait 120 \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args "-B $WD,$GENOME,$BISMARK_IDX,$TMP_FOLDER" \
    --jobscript jobscripts/narvi.sh \
    --rerun-incomplete
