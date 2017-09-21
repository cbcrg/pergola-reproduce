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
### Gviz plot of raw jaaba variables derived from fly trajectory         ###
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
      gviz_pergola_bedgraph
      
      Arguments:
      --path2variables=path_to_bedgr_files - character
      --variable_name=someValue            - character, variable name
      --image_format=\"image_format\"
      --help                               - print this text
      
      Example:
      ./melanogaster_gviz_var.R --path_bed_files=\"path_bed_files\" --image_format=\"image_format\" \n")      
      
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

# plot image format
{
    if (is.null (argsL$image_format))
    {
        image_format <- "tiff"
        warning ("[Warning]: format for plots not provided, default tiff")
    }
    else
    {
        image_format <- argsL$image_format
    }
}

## Loading libraries
# library(utils)
# source("https://bioconductor.org/biocLite.R")
# library(gtools) #mixedsort

## bedGraph files to GRanges
library("GenomicRanges")
library("rtracklayer")
library("Gviz")

##############################
## bedgraph row measures files
bedgraph.var_files <- list.files(path=path2bedg_files, pattern="^tr*", full.names=TRUE)
# bedgraph.var_files

bedgraph_tracks <- lapply(bedgraph.var_files, function (bedgraph) {
#   id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
  id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bedgraph)
#   id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedgraph)))
#   print (id)
  bedgraph_GR <- import(bedgraph, format = "bedGraph")
  tr <- DataTrack(bedgraph_GR, type="a", ylim = c(0, 50),
            #background.title = l_gr_color[[i_group_exp]],           
            showAxis = F, name = id)
  return (tr)
} )

names(bedgraph_tracks) <- as.numeric(gsub(".+tr_(\\d+)(_.+$)", "\\1", bedgraph.var_files))
bedgraph_tracks <- bedgraph_tracks[as.character(1:length(bedgraph.var_files))]

plot_name <- "gviz_jaaba_var"

{
    if (image_format == 'tiff' | image_format == 'tif') {
        tiff(paste(plot_name, ".", image_format, sep=""),  width = 45 , height = 34, units = "cm", res=300)
    }
    else if (image_format == 'pdf') {        
        pdf(paste(plot_name, ".", image_format, sep=""),  width = 45 , height = 34)
    }
    else if (image_format == 'png') {        
        png(paste(plot_name, ".", image_format, sep=""), width = 45 , height = 34, units = "cm", res=300)
    }
    else {
        stop (paste("Unknow image file format:", image_format, sep=" "))
    }
}

plotTracks(bedgraph_tracks, stacking="dense", collapse=FALSE, shape = "box", col="red")

dev.off()
