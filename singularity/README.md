[![hackmd-github-sync-badge](https://hackmd.io/q4nlB0nJRLuvu7x96Gri2A/badge)](https://hackmd.io/q4nlB0nJRLuvu7x96Gri2A)

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

A pre-compiled container image can be pulled from the Singularity Library. This is the recommended way to use this container:

```bash
 singularity pull library://ftabaro/default/methylsnake 
 ```

## Compiling the container 

The container images can be built locally with:

```bash
sudo singularity build base.sif Singularity.base
sudo singularity build methylkit.sif Singularity.methylkit
```

The `Singularity.methylkit` recipe expects a container called `base.sif` in the same directory as the recipe.

###### tags: `singularity` `methylsnake` `DNA-methylation` `rrbs`
