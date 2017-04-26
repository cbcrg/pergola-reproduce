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
        --path2variables=someValue - character, path to read bedGraph files
        --path2scores=someValue    - character, path to read bed files        
        --variable_name=someValue  - character, variable name
        --help                     - print this text
        
        Example:
        ./sushi_plot_jaaba.R --path2variables=\"/foo/variables\" --path2scores=\"/foo/scores\" --variable_name=\"var_name\"\n")
    
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

## Loading libraries
library(utils)
source("https://bioconductor.org/biocLite.R")
library('Sushi')
library(gtools) #mixedsort

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

##############################
## bedgraph row measures files
base_folder <- path2bedg_files
chase.bedgraph.variable.files <- mixedsort(list.files(base_folder, pattern=paste("tr.*", variable_name, "*.bedGraph$", sep=""), full.names=TRUE))

data_bedgraph_variable <- lapply(chase.bedgraph.variable.files, function (bedg) { 
    name_id <- gsub("tr_", "", gsub(paste("_dt_", variable_name, ".bedGraph", sep=""), "", basename(bedg)))
    bedg_tbl <- read.csv(file=bedg, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    bedg_tbl$name <- as.numeric(name_id)
    return (bedg_tbl)
})

chrom            = "chr1"
chromstart       = 0
chromend         = 25000

title <- paste("  ", variable_name, sep="")

png(paste("sushi_jaaba_scores_annot_", variable_name, ".png", sep=""))

split.screen (c(2, 1)) 

## adding a n empty plots for title
n=3
split.screen(c(length(data_bedgraph_variable)+n, 1), screen = 1)

screen(1)

par(mar=c(0.1,1,2,0.1))

plot(1, type="n", axes=F, xlab="", ylab="")

labelplot("A ", title, letteradj=-.025)

i=3+n
j=1

for (bedg in data_bedgraph_variable) {
    screen( i )
    par(mar=c(0.1,1,0.1,0.1))
    plotBedgraph(bedg, chrom, chromstart, chromend, transparency=.50,
                 color=SushiColors(2)(2)[1])
    axis(side=2, las=2, tcl=.2)    
    mtext(j, side=2, cex=0.5)
    i=i+1
    j=j+1
}

screen( 2 )
par(mar=c(0.1,1,2,0.1))
plotBed (beddata = data_bed.df, chrom = chrom,
         chromstart = chromstart, chromend = chromend,
         rownumber  = data_bed.df$row, type = "region",
         color = data_bed.df$color, row ="given",
         plotbg ="grey95", rowlabels = rev(unique(data_bed.df$name)),
         rowlabelcol = unique(data_bed.df$color), rowlabelcex = 0.75)
#rowlabelcex = 0.75)
labelplot("B ","  Annotations", letteradj=-.025)

dev.off()