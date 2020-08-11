[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/4638) [![hackmd-github-sync-badge](https://hackmd.io/q4nlB0nJRLuvu7x96Gri2A/badge)](https://hackmd.io/q4nlB0nJRLuvu7x96Gri2A)

# MethylSnake Singularity container

This is the recipe file for the Singularity container for the MethylSnake pipeline. It contains all the tools required to run it:
- samtools
- bowtie2
- TrimGalore!
- Bismark
- R 4.0.2 + Bioconductor + MethylKit
- gtfToGenePred + genePredToBed

## Pulling the container

A pre-compiled container image can be pulled from SingularityHub. This is the recommended way to use this container:

```bash
singularity pull shub://ftabaro/MethylSnake:latest
```

## Compiling the container 

The container image can be built locally with:

```bash
singularity build MethylSnake.sif Singularity
```

###### tags: `singularity` `methylsnake` `DNA-methylation` `rrbs`
