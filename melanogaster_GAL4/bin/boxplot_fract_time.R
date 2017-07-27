#!/usr/bin/env Rscript

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
############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Aug 2017                        ###
############################################################################
### Boxplot comparing fraction of time performing a behavior of several  ###
### drosophila strains                                                   ###
############################################################################

#####################
## VARIABLES
## Reading arguments
args <- commandArgs (TRUE) #if not it doesn't start to count correctly

## Default setting when no arguments passed
if ( length(args) < 1) {
    args <- c("--help")
}

## Help section
if("--help" %in% args) {
    cat("
        boxplot_fract_time

        Arguments:
        --path2tbl_fr_time=someValue - character, path to read tbl files
        --help                     - print this text

        Example:
        ./sushi_plot_jaaba.R --path2tbl_fr_time=\"/foo/tbl.txt\"\n")

    q (save="no")
}

## Use to parse arguments beginning by --
parseArgs <- function(x)
{
    strsplit (sub ("^--", "", x), "=")
}

## Parsing arguments
argsDF <- as.data.frame (do.call("rbind", parseArgs(args)))
argsL <- as.list (as.character(argsDF$V2))
names (argsL) <- argsDF$V1

# path to variables bedgraph files
{
    if (is.null (argsL$path2tbl_fr_time))
    {
        stop ("[FATAL]: Path to files variables bedgraph files is mandatory")
    }
    else
    {
        path2tbl_fr_time <- argsL$path2tbl_fr_time
    }
}

## Loading libraries
library (ggplot2)

## Loading parameters for plots
source("https://gist.githubusercontent.com/JoseEspinosa/307430167b05e408ac071b8724701abf/raw/06b26f764953ceea53d334190d7b736308e5141d/param_plot_publication.R")

tbl_frac_time <- read.csv(path2tbl_fr_time, sep="\t", header=F)
colnames(tbl_frac_time) <- c("fraction", "group")


## colors
cbb_palette <- c("#D55E00", "#0072B2", "#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
variable <- "boxplot_fract_time"
behavior <- gsub ("_.*$", "", tbl_frac_time[1,2])
name_out <- paste (variable, "_", behavior, ".", "pdf", sep="")
tbl_frac_time$strain <- gsub("^.*_", "", as.character(tbl_frac_time$group))

ggplot(tbl_frac_time, aes(strain, fraction, fill=strain)) + geom_boxplot(notch=FALSE) +
  labs (#title = "Jaaba annotated vs. non-annotated intervals\n",
    y = paste("Fraction time ", behavior, "\n", sep=""), x = "\n") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = cbb_palette) +
  theme(legend.position="none") +
  ylim(0,max(tbl_frac_time$fraction))

ggsave (file=name_out)
