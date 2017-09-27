#  Copyright (c) 2014-2017, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2017, Jose Espinosa-Carrasco and the respective authors.
#
#  This file is part of Pergola.
#
#  Pergola is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Pergola is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################
# Dockerfile to build Pergola and all dependencies to reproduce 
# pergola paper results
################################################################

FROM r-base:3.3.2

MAINTAINER Jose Espinosa-Carrasco <espinosacarrascoj@gmail.com>

RUN apt-get update && \
    apt-get install --fix-missing -y \
    sudo \
    python \
    python-pip \
    bedtools

RUN pip install pybedtools

## pergola installation
COPY pergola/pergola /pergola/pergola
COPY pergola/requirements.txt /pergola/
COPY pergola/setup.py /pergola/
COPY pergola/README.md /pergola/

RUN pip install -r /pergola/requirements.txt && \
    pip install cython && \
    pip install h5py && \
    apt-get install -y python-scipy && \
    cd pergola && python setup.py install

## install R dependencies
RUN apt-get update && \
    apt-get install --fix-missing -y \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
    libxml2-dev \
    python-dev \
    python-distribute \
    gfortran \
    bedtools \
    libhdf5-dev

## Install R packages
RUN R -e "install.packages(c('Sushi', 'shiny', 'rmarkdown', 'ggplot2', 'XML', 'Rcurl','cowplot', 'dplyr', 'survival', 'gridExtra', 'utils', 'gutils', 'gtools', 'ggrepel', 'extrafont', 'devtools'), repos='http://cran.rstudio.com/')" \
&& Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges"); biocLite("rtracklayer"); biocLite("Sushi");'

# version of Gviz modified to show time instead of genomics units
RUN R -e  'devtools::install_github("JoseEspinosa/Gviz")'

## version of Gviz modified to show fps instead of genomics units
# RUN R -e  'devtools::install_github("JoseEspinosa/Gviz", ref = "fps")'

