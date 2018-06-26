#  Copyright (c) 2014-2018, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2018, Jose Espinosa-Carrasco and the respective authors.
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
    python-dev \
    python-distribute \
    python-pip \
    bedtools

## install R dependencies
RUN apt-get update && \
	apt-get install --fix-missing -y \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libxt-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev

## Install R packages
RUN R -e "install.packages(c('Sushi', 'shiny', 'rmarkdown', 'ggplot2', 'XML', 'Rcurl','cowplot', 'dplyr', 'survival', 'gridExtra', 'utils', 'gutils', 'gtools', 'ggrepel', 'extrafont', 'devtools', 'withr'), repos='http://cran.rstudio.com/')" \
&& Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges"); biocLite("rtracklayer"); biocLite("Sushi");'

# Two versions of gviz, display time axis as fps or time units
RUN bash -c 'mkdir -p /gviz/fps; mkdir -p /gviz/time'

# version of Gviz modified to show time units instead of genomics units
RUN R -e  'withr::with_libpaths(new = "/gviz/time", devtools::install_github("JoseEspinosa/Gviz"))'

## version of Gviz modified to show fps instead of genomics units
RUN R -e  'withr::with_libpaths(new = "/gviz/fps", devtools::install_github("JoseEspinosa/Gviz", ref = "fps"))'

RUN pip install tables

COPY pergola/requirements.txt /pergola/
RUN pip install -r /pergola/requirements.txt
RUN pip install cython

# I need this to avoid the broken package list apt
RUN rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y --no-install-recommends git \
                                               libssl-dev \
                                               openssl \
                                               mysql-client-5.7 \
                                               mysql-client-core-5.7 \
                                               libmysqlclient-dev

# Compile and install kentUtils
RUN cd /tmp && \
    git clone https://github.com/ENCODE-DCC/kentUtils.git && \
    cd kentUtils && \
    git checkout v302.1.0 && \
    make && \
    cp -rp bin/* /usr/local/bin && \
    cd .. && rm -rf kentUtils

# Compile and install bwtool
RUN git clone https://github.com/CRG-Barcelona/libbeato.git && \
    git clone https://github.com/CRG-Barcelona/bwtool.git && \
    cd libbeato/  &&  ./configure && \
    make && \
    sudo make install && \
    cd ../bwtool/  && ./configure && \
    make && \
    sudo make install

# Install deeptools
# RUN pip install deeptools==2.5.7
RUN pip install deeptools==3.0.2

# Install java jdk and chromHMM 1.15
RUN sudo apt-get install --fix-missing -y default-jdk
RUN git clone --branch v1.15 https://github.com/jernst98/ChromHMM --single-branch

# Update path with the scripts
ENV PATH="/ChromHMM:${PATH}"

# Install pergola
COPY pergola/pergola /pergola/pergola
COPY pergola/requirements.txt /pergola/
COPY pergola/setup.py /pergola/
COPY pergola/README.md /pergola/

RUN cd pergola && python setup.py install

# Reinstall matplotlib otherwise deeptools crashes
RUN pip install matplotlib==2.0.2
