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
### Jose Espinosa-Carrasco CB-CRG Group. Dec 2016                        ###
############################################################################
### Sushi plot of raw jaaba variables derived from the trajectory        ###
###                                                                      ### 
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
      sushi_plot_jaaba
      
      Arguments:
      --path2variables=someValue - character, path to read bedGraph files    
      --variable_name=someValue  - character, variable name
      --help                     - print this text
      
      Example:
      ./sushi_plot_jaaba.R --path2variables=\"/foo/variables\" --variable_name=\"var_name\"\n")
  
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
  if (is.null (argsL$path2variables)) 
  {
    stop ("[FATAL]: Path to files variables bedgraph files is mandatory")
  }
  else
  {
    path2bedg_files <- argsL$path2variables
  }
}

# variable name
{
  if (is.null (argsL$variable_name)) 
  {
    print ("[WARNING]: Variable to plot is set to \"default=Speed\"\n")
    variable_name <- "velmag"
  }
  else
  {        
    variable_name <- argsL$variable_name
  }
}

# behavior_group
{
  if (is.null (argsL$behavior_strain))
  {
    print ("[WARNING]: behavior_group is set to \"default=behavior_strain\"\n")
    behavior_strain <- "behavior_strain"
  }
  else
  {
    behavior_strain <- argsL$behavior_strain
  }
}

## Loading libraries
library(utils)
source("https://bioconductor.org/biocLite.R")
library('Sushi')
library(gtools) #mixedsort

##############################
## bedgraph row measures files
# path2bedg_files <- "/Users/jespinosa/git/pergola-paper-reproduce/melanogaster_GAL4/results/results_bedGr/"
# variable_name <- "velmag"
base_folder <- path2bedg_files
chase.bedgraph.variable.files <- mixedsort(list.files(base_folder, pattern=paste("values.*", variable_name, ".bedGraph$", sep=""), full.names=TRUE))
chase.bedgraph.variable.files.comp <- mixedsort(list.files(base_folder, pattern=paste("values.*", variable_name, "*.comp.*.bedGraph$", sep=""), full.names=TRUE))

# min_bedGraph <- 0
# max_bedGraph <- -100000
# 
# round_up <- function(x,to=10) {
#   to*(x%/%to + as.logical(x%%to))
# }

data_bedgraph_variable <- lapply(chase.bedgraph.variable.files, function (bedg) { 
  name_id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedg)))  
  bedg_tbl <- read.csv(file=bedg, header=FALSE, sep="\t", stringsAsFactors=FALSE)
#   min_bedGraph <<- round(floor(min (min_bedGraph, min(bedg_tbl$V4))),-1)
#   max_bedGraph <<- round_up(max (max_bedGraph, max(bedg_tbl$V4)), 5)
  bedg_tbl$name <- as.numeric(name_id)
  return (bedg_tbl)
})

data_bedgraph_variable_comp <- lapply(chase.bedgraph.variable.files.comp, function (bedg) { 
  name_id <- gsub("values_", "", gsub(paste("_", variable_name, ".comp.bedGraph", sep=""), "", basename(bedg)))
  bedg_tbl <- read.csv(file=bedg, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  bedg_tbl$name <- as.numeric(name_id)
#   min_bedGraph <<- round(floor(min (min_bedGraph, min(bedg_tbl$V4))),-1)
#   max_bedGraph <<- round_up(max (max_bedGraph, max(bedg_tbl$V4)), 5)
  return (bedg_tbl)
})

ylim = boxplot.stats(do.call(rbind, data_bedgraph_variable)$V4)$stats[c(1, 5)]
ylim_comp = boxplot.stats(do.call(rbind, data_bedgraph_variable_comp)$V4)$stats[c(1, 5)]
ranges_plot <- c(min(ylim, ylim_comp), max(ylim, ylim_comp))

chrom            = "chr1"
chromstart       = 0
chromend         = 25000

# title <- paste("  ", variable_name, sep="")

# png(paste("sushi_highlight", variable_name, "_", behavior_strain, ".png", sep=""))
pdf ( paste("sushi_highlight", variable_name, "_", behavior_strain, ".pdf", sep="") , height=10, width=20)

# split.screen (c(2, 1)) 

## adding a n empty plots for title
n=3
# split.screen(c(length(data_bedgraph_variable)+n, 1), screen = 1)
split.screen(c(length(data_bedgraph_variable), 1))
screen(1)

par(mar=c(0.1,1,2,0.1))

# plot(1, type="n", axes=F, xlab="", ylab="")

# labelplot("A ", title, letteradj=-.025)

# i=3+n
i=1
j=1

for (bedg_i in seq_along(data_bedgraph_variable)) {
  screen( i )
  #   par(mar=c(0.1,1,0.1,0.1))
  par(mar=c(0.1, 2, 0.1, 2))
  plotBedgraph(data_bedgraph_variable_comp[[bedg_i]], chrom, chromstart, chromend,
               #                transparency=.10, 
               transparency=.50,
               color=SushiColors(2)(2)[2], range=ranges_plot)
  plotBedgraph(data_bedgraph_variable[[bedg_i]], chrom, chromstart, chromend, transparency=.50,
               color=SushiColors(2)(2)[1], overlay=TRUE)# , rescaleoverlay=TRUE)  
  #   axis(side=2,las=2,tcl=.2)
  axis(side=2, las=1, tcl=.2, cex.axis=0.5) 
  #   mtext("Speed (mm/s)\n",side=2,line=1.75,cex=1,font=2)
  #   mtext(j, side=2, cex=0.5)
  mtext(j, side=2, cex=0.5)
  i=i+1
  j=j+1
}

dev.off()
