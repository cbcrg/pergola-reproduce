# melanogaster_GAL4-Pergola-Reproduce.nf

This repository contains the software, scripts and data to reproduce the results corresponding to the *D. melanogaster* GAL4 line experiment of the Pergola paper.

If you have not install yet [docker](https://www.docker.com/) and [nextflow](https://www.nextflow.io/), follow this [intructions](../README.md)

## Clone the repository

```bash
git clone --recursive https://github.com/JoseEspinosa/pergola-paper-reproduce.git
cd pergola-paper-reproduce/melanogaster_GAL4
```

**Note**: If you have previously download the repository, then you only need to go to the ``melanogaster_GAL4`` folder.

## Data

Data is publicly available in [Zenodo](https://zenodo.org/) as a compressed tarball.

Data can be downloaded and uncompressed using the following command:

```bash
mkdir data
wget -O- https://zenodo.org/record/582475/files/melanogaster_GAL4.tar.gz | tar xz -C data
```

#### Original Data Sources
The original data can be downloaded from the following [source](http://sourceforge.net/projects/jaaba/files/Sample%20Data/sampledata_v0.5.zip/download).

## Pull docker image
Pull the Docker image use for processing data with Pergola (Pergola and its dependencies installed)

```bash
docker pull pergola/pergola@sha256:f7208e45e761dc0cfd3e3915237eb1a96eead6dfa9c8f3a5b2414de9b8df3a3d
```

## Run nextflow pipeline
Once data is downloaded, it is possible to reproduce all the results using this command:

```bash
NXF_VER=0.26.1 nextflow run melanogaster_GAL4-Pergola-Reproduce.nf \ 
    --scores='data/scores/scores_chase_*.mat' \
    --var_dir='data/perframe_*' \
    --variables="velmag" \
    --mappings='data/mappings/jaaba2pergola.txt' \
    -with-docker
```

## Sushi visualization

The previous command generates a results folder that contains the plots used in the paper figure:

* A Sushi plot the intervals corresponding to chasing behavior annotated using [Jaaba](http://jaaba.sourceforge.net/).
* A Sushi plot of each of the variables derived from the flies trajectories using [CTRAX](http://ctrax.sourceforge.net/). This variables are used by JAABA to train the classifier and predict the intervals annotated as chasing behavior. 
* A volcano plot that shows the fold change of the variables from intervals annotated as chansing behavior and the rest of the trajectory, the magnitude of the change, and how significant it is.
* A box plot for each of the variables comparing chase annotated intervals and non-annotated.
* Besides in the results folder you can find all the intermediate files used for the rendering of the data.

**Note**: Sushi is a R/Bioconductor package, designed for the creation of publication-quality plots for genomic visualizations, more info [here](https://www.bioconductor.org/packages/release/bioc/html/Sushi.html).
