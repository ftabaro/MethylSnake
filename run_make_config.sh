python scripts/make_config.py --config-path /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/config.yaml \
  --wd /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma \
  --genome-path /bmt-data/genomics/reference/hg38.fa \
  --bismark-index-path /bmt-data/genomics/reference/bismark/hg38 \
  --annotation-file /bmt-data/genomics/reference/Ensembl/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz \
  --sample-sheet /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/sample_sheet.csv \
  --tmp-folder /scratch/ft413468/rrbs-sharma \
  --log-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/log \
  --reads-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/reads \
  --trimmed-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/trimmed \
  --alignments-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/alignments \
  --reports-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/reports \
  --nucleotide-stats-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/nucl_stat \
  --rdata-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/Rdata \
  --pictures-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/pictures \
  --tables-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/tables \
  --mate1-pattern _R1_001 \
  --mate2-pattern _R2_001 \
  --fastq-extension .fastq.gz 
#  --environments-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/rrbs-pipeline/environments 
#  --methylkitdb-folder /lustre/bmt-data/genomics/projects/RRBSdata_190047_Sharma/methylkitdb \

