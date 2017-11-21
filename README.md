# pergola-paper-reproduce

This repository contains the software, scripts and data to reproduce the results and figures decribed in the Pergola publication.

* To reproduce *C. elegans* results go to [worm](celegans_n2_unc16/README.md) 

* To reproduce *D. melanogaster* results go to [fly](melanogaster_GAL4/README.md)

* To reproduce *M. musculus* *(mice) results go to [mouse](cb1_mice/README.md) 

In order to avoid any possible problem reproducing the results we recommend to install [docker](https://www.docker.com/) following the instructions on this [link](https://docs.docker.com/engine/installation/). Using the images we built guarantees that all the processes runs exactly in the same system environment and using the same software and libraries versions, you can read more about it [here](https://peerj.com/articles/1273/).

The pipelines followed in each of the three cases were coded using [nextflow](https://www.nextflow.io/). Nextflow scripts allow to combine several scripting languages and the use of docker container out of the box. To intall nextflow just execute the following command. 

```
curl -fsSL get.nextflow.io | bash
```

