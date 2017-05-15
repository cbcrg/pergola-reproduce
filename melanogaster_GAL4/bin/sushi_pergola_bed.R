#!/usr/bin/env Rscript

#  Copyright (c) 2014-2016, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2016, Jose Espinosa-Carrasco and the respective authors.
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
      --path2scores=someValue    - character, path to read bed files              
      --help                     - print this text
      
      Example:
      ./sushi_pergola_bed.R --path2scores=\"/foo/scores\" \n")
  
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

# path to score bed files
{
  if (is.null (argsL$path2scores)) 
  {
    stop ("[FATAL]: Path to files score bed files is mandatory")
  }
  else
  {
    path2bed_files <- argsL$path2scores
  }
}

## Loading libraries
library(utils)
source("https://bioconductor.org/biocLite.R")
library('Sushi')
library(gtools) #mixedsort

# path2bed_files <- "/Users/jespinosa/git/pergola-paper-reproduce/melanogaster_GAL4/results/results_score/"
base_folder <- path2bed_files
chase.bed.files <- mixedsort(list.files(base_folder, pattern="tr.*.bed$", full.names=TRUE))
data_bed <- lapply(chase.bed.files, function (bed, dir=direction) { 
  name_id <- gsub("tr_", "", gsub("_dt_chase.bed", "", basename(bed)))
  bed_tbl <- read.csv(file=bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  bed_tbl$name <- as.numeric(name_id)
  return (bed_tbl)
})

## ROW are reversed!!! otherwise first row is 20
## also rowlabels must be reversed for the correct labelling in plot
data_bed.df <- do.call(rbind, data_bed)
data_bed.df <- data_bed.df[order(data_bed.df$name),] 
data_bed.df <- transform(data_bed.df, row=match(data_bed.df$name, rev(unique(data_bed.df$name))))
data_bed.df$color <- sapply(strsplit(as.character(data_bed.df$V9), ","), function(x) {
  rgb(x[1], x[2], x[3], m=255)
})   

chrom            = "chr1"
chromstart       = 0
chromend         = 25000

# png(paste("sushi_jaaba_annot", ".png", sep=""))
pdf ( paste("sushi_jaaba_annot", ".pdf", sep="") , height=10, width=20)

# par(mar=c(0.1,1,2,0.1))
par(mar=c(1, 2, 1, 2))

plotBed (beddata = data_bed.df, chrom = chrom,
         chromstart = chromstart, chromend = chromend,
         rownumber  = data_bed.df$row, type = "region",
         color = data_bed.df$color, row ="given",
         plotbg ="grey95", rowlabels = rev(unique(data_bed.df$name)),
#          rowlabelcol = unique(data_bed.df$color), rowlabelcex = 1)
         rowlabelcol = "black", rowlabelcex = 1)
#rowlabelcex = 0.75)
# labelplot("B ","  Chase", letteradj=-.025)

dev.off()