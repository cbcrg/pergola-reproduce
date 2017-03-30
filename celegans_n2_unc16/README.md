# N2_vs_unc16_motions-Pergola-Reproduce.nf

This repository contains the software, scripts and data to reproduce the results corresponding to the C.elegans data of the Pergola paper.

If you have not install yet [docker](https://www.docker.com/) and [nextflow](https://www.nextflow.io/), follow this [intructions](../README.md)

## Pull docker image
Pull the Docker image use for processing data with Pergola (Pergola and its dependencies installed)

```
docker pull pergola/pergola@sha256:6d032d23fd90317ee3ac564497b1fd7f204ffc641ee009b937846fe7c959834f
```

## Clone the repository

```
git clone --recursive https://github.com/JoseEspinosa/pergola-paper-reproduce.git
cd pergola-paper-reproduce/celegans_n2_unc16
```

## Data
Data is publicly available in [Zenodo](https://zenodo.org/) as a compressed tarball.

Data can be downloaded and uncompressed using the following command:

```
mkdir data
wget -O- https://zenodo.org/record/400948/files/celegans_unc16_N2.tar.gz | tar xz -C data # todo modify this part 
```

## Run nextflow pipeline
Once data is downloaded, it is possible to reproduce all the results using this command:

```
NXF_VER=0.24.1 ./N2_vs_unc16_motions-Pergola-Reproduce.nf --strain1_trackings 'data/unc_16/*.mat' --strain2_trackings 'data/N2/*.mat' \
  --mappings_speed 'data/mappings/worms_speed2p.txt' \
	--mappings_bed 'data/mappings/bed2pergola.txt' \
	--mappings_motion data/mappings/worms_motion2p.txt -with-docker
```	


## Online visualization

TODO