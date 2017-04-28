# melanogaster_GAL4-Pergola-Reproduce.nf

This repository contains the software, scripts and data to reproduce the results corresponding to the *D. melanogaster* GAL4 line experiment of the Pergola paper.

If you have not install yet [docker](https://www.docker.com/) and [nextflow](https://www.nextflow.io/), follow this [intructions](../README.md)

## Data processing

### Pull docker image
Pull the Docker image use for processing data with Pergola (Pergola and its dependencies installed)

```
docker pull pergola/pergola@
```








### Clone the repository

```
git clone --recursive https://github.com/JoseEspinosa/pergola-paper-reproduce.git
cd pergola-paper-reproduce/melanogaster_GAL4
```

### Data

Data is publicly available in [Zenodo](https://zenodo.org/) as a compressed tarball.

Data can be downloaded and uncompressed using the following command:

```
mkdir data
# TODO wget -O- https://zenodo.org/record/398779/files/CB1_mice.tar.gz | tar xz -C data
wget -O-  https://zenodo.org/record/398779/files/melanogaster_gal4.tar.gz | tar xz -C data
```

### Run nextflow pipeline
Once data is downloaded, it is possible to reproduce all the results using this command:

```
NXF_VER=0.24.0 ./melanogaster_GAL4-Pergola-Reproduce.nf --scores='data/scores/*.mat' \
                                                        --var_dir='data/perframe/' \
                                                        --variables="velmag" \
                                                        --mappings='data/jaaba2pergola.txt' \
                                                        -with-docker
```


## Sushi visualization


The previous command generates a results folder that contains the figures named sushi_jaaba_scores_annot_velmag.png. The plot
shows the 'velmag' (speed of the center of rotation) variable derived from the trajectory of the flies and the intervals
corresponding to chasing behavior annotated using [Jaaba](http://jaaba.sourceforge.net/). This figure is genereated using
Sushi an R/Bioconductor package, a package designed for the creation of publication-quality plots for genomic visualizations.
Besides in the results folder you can find all the intermediate files used for the rendering of the data.


