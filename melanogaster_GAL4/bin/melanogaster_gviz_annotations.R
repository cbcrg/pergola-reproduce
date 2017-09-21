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
#############################################################################
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. Aug 2017                   ###
#############################################################################
### Mice data set visualization using Gviz                                ###
#############################################################################

##Getting HOME directory
home <- Sys.getenv("HOME")

### Execution example
## Rscript melanogaster_gviz_visualization.R --path_bed_files="path_to_bed_files"
library("ggplot2")
library("Gviz")

## bed files to GRanges
library("GenomicRanges")
library("rtracklayer")

# Loading params plot:
source("https://raw.githubusercontent.com/cbcrg/mwm/master/lib/R/plot_param_public.R")

#####################
### VARIABLES
#Reading arguments
args <- commandArgs (TRUE) #if not it doesn't start to count correctly

## Default setting when no arguments passed
if ( length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      density_heatmaps_tracking.R
      
      Arguments:
      --path_bed_files=path_to_bed_files - character
      --image_format=image_format        - character
      --help                             - print this text
      
      Example:      
      ./melanogaster_gviz_annotations.R --path_bed_files=\"path_bed_files\" --image_format=\"image_format\" \n")
  
  q (save="no")
}

# Use to parse arguments beginning by --
parseArgs <- function(x)
{
  strsplit (sub ("^--", "", x), "=")
}

#Parsing arguments
argsDF <- as.data.frame (do.call("rbind", parseArgs(args)))
argsL <- as.list (as.character(argsDF$V2))
names (argsL) <- argsDF$V1

# All arguments are mandatory
{
  if (is.null (argsL$path_bed_files))
  {
    stop ("[FATAL]: path to bed files not provided", stderr())
  }
  else
  {
    path_bed_files <- argsL$path_bed_files
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

#############################
## Read files bed files
bed_files <- list.files(path=path_bed_files, pattern="^tr*", full.names=TRUE)

bed_tracks <- lapply(bed_files, function (bed) {
                id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
                bed_GR <- import(bed, format = "bed")
                tr <- AnnotationTrack(bed_GR, name = paste ("", id, sep=""),
                      fill=bed_GR$itemRgb)
#                       background.title = c)  
  names(tr) <- id
  return (tr)
  } )

names(bed_tracks) <- as.numeric(gsub(".+tr_(\\d+)(_.+$)", "\\1", bed_files))
bed_tracks <- bed_tracks[as.character(1:length(bed_files))]

## Plot
plot_name <- "gviz_jaaba_annot"

{
    if (image_format == 'tiff' | image_format == 'tif') {
        tiff(paste(plot_name, ".", image_format, sep=""), width = 45 , height = 34, units = "cm", res=300)
    }
    else if (image_format == 'pdf') {        
        pdf(paste(plot_name, ".", image_format, sep=""), height=45, width=34)        
    }
    else if (image_format == 'png') {        
        png(paste(plot_name, ".", image_format, sep=""),  width = 45 , height = 34, units = "cm", res=300) 
    }
    else {
        stop (paste("Unknow image file format:", image_format, sep=" "))
    }
}

# tiff(name_file, width = 45 , height = 34, units = "cm", res=300)
plotTracks(bed_tracks, stacking="dense", from=0, collapse=FALSE, shape = "box", col=NULL) 

dev.off()
