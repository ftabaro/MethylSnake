[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4638) [![hackmd-github-sync-badge](https://hackmd.io/q4nlB0nJRLuvu7x96Gri2A/badge)](https://hackmd.io/q4nlB0nJRLuvu7x96Gri2A)

# MethylSnake Singularity container

This folder holds two recipes to run the MethylSnake pipeline. The `Singularity.base` recipe creates a container with basic tools (Bowtie2, TrimGalore!, Bismark, etc.), while the `Singularity.methylkit` recipe installs R and the required packages on top of the base container. To run the pipeline only the `methylkit` container is required. 

The final container (`methylkit`) contain:
- samtools
- bowtie2
- TrimGalore!
- Bismark
- R 4.0.2 + Bioconductor + MethylKit
- gtfToGenePred + genePredToBed

## Pulling the container

A pre-compiled container image can be pulled from SingularityHub. This is the recommended way to use this container:

```bash
singularity pull shub://ftabaro/MethylSnake:methylkit
```

## Compiling the container 

The container images can be built locally with:

```bash
sudo singularity build base.sif Singularity.base
sudo singularity build methylkit.sif Singularity.methylkit
```

###### tags: `singularity` `methylsnake` `DNA-methylation` `rrbs`
