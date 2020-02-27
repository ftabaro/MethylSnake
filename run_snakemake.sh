#!/bin/bash
if [ $# -eq 0 ]; then
    echo "Config file is required."
    exit 1
fi

LOG_FOLDER=$(grep log_folder "$1" | awk -F":" '{print $2}')
mkdir -p "$LOG_FOLDER"

WD=$(grep wd "$1" | awk -F":" '{print $2}')
GENOME=$(grep genome_path "$1" | awk -F":" '{print $2}' |xargs dirname)
BISMARK_IDX=$(grep bismark_index_path "$1" | awk -F":" '{print $2}' |xargs dirname)
ANNOTATION=$(grep annotation_file "$1" | awk -F":" '{print $2}' |xargs dirname)
TMP_FOLDER=$(grep tmp_folder "$1" | awk -F":" '{print $2}')

PATHS="$WD,$ANNOTATION,$GENOME,$BISMARK_IDX,$TMP_FOLDER"
PATHS=$(echo $PATHS |sed 's/ //g')


SBATCH="sbatch -J {cluster.jobname} --partition={cluster.partition} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus} -e ${LOG_FOLDER%/}/{cluster.log} -o ${LOG_FOLDER%/}/{cluster.log}"

snakemake -r -p -j 100 --configfile "$1" \
    --cluster "$SBATCH" \
    --cluster-config cluster-config/cluster.json \
    --latency-wait 120 \
    --rerun-incomplete \
    --use-singularity \
    --singularity-args "--bind $PATHS" \
    --rerun-incomplete
#    --jobscript jobscripts/narvi.sh \
