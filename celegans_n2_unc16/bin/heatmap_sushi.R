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
### Heatmap Sushi plot of C.elegans mid body speed                       ###
###                                                                      ### 
############################################################################

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
        --path_bedgr=path_to_files  - character  
        --image_format=image_format - character
        --help                      - print this text
        
        Example:
        ./plot_speed_motion_mean.R --path_bedgr=\"path_to_files\" --image_format=\"image_format\" \n")
    
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
    if (is.null (argsL$path_bedgr)) 
    {
        stop ("[FATAL]: path_bedgr arg is mandatory")
    }
    else
    {
        path2bedg_files <- argsL$path_bedgr
    }
}

{
    if (is.null (argsL$image_format))
    {
        image_format <- ".tiff"
        warning ("[Warning]: format for plots not provided, default tiff")        
    }
    else
    {
        image_format <- argsL$image_format
    }
}

# ## Loading libraries
library('Sushi')

## bedgraph behavioral measures files
# path2bedg_files <- "/Users/jespinosa/git/pergola-paper-reproduce/celegans_n2_unc16/results/results_bedgr1/"
# path2bedg_files<- "/Users/jespinosa/scratch/b4/a5db1241fd345a668a8e8234ee595b/results_bedgr1"
#############################
## Read bedGraph files
bedg_files <- list.files(path=path2bedg_files, pattern="^bedg*", full.names=TRUE)

df_order <- as.data.frame(bedg_files)
df_order$id_num <- as.numeric(sub("bedg.str[1-2].(.*?)\\.bedGraph", "\\1", basename(bedg_files)))
bedg_files <- as.vector(df_order[with(df_order, order(id_num)), "bedg_files"])

unite_scale <- function (v, min=-400, max=400) {
    if (v < -400 & v != -10000) { return (1) }
    else if (v > 400 ) { return (1) }
    else if (v == -10000) { return (0) } 
    else if (v < 0) { return ((v - min) / (0 - min)) } 
    return ((v - 0) / (max - 0))
}

data_bedgraph_variable <- lapply(bedg_files, function (bedg) { 
    
    name_id <- gsub("bedg.str1.", "", gsub(".bedGraph", "", basename(bedg)))  
    bedg_tbl <- read.csv(file=bedg, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    bedg_tbl$name <- as.numeric(name_id)    
    bedg_tbl$color <- ifelse(bedg_tbl$V4 > 0, opaque("red", transparency=sapply(bedg_tbl$V4, unite_scale)), opaque("blue", transparency=sapply(bedg_tbl$V4, unite_scale)))
    
    return (bedg_tbl)
})

chrom            = "chr1"
chromstart       = 0
chromend         = 29002

{
    if (image_format == 'tiff' | image_format == 'tif') {
        tiff(paste("sushi_var", ".", image_format, sep="") , height=10, width=20, units="cm", res=300)
        size_lab <- 0.3
    }
    else if (image_format == 'pdf') {        
        pdf(paste("sushi_var", ".", image_format, sep="") , height=10, width=20)
        size_lab <- 0.5
    }
    else if (image_format == 'png') {        
        png(paste("sushi_var", ".", image_format, sep=""))
        size_lab <- 0.3
    }
    else {
        stop (paste("Unknow image file format:", image_format, sep=" "))
    }
}

## adding a n empty plots for title
n=3
# split.screen(c(length(data_bedgraph_variable)+n, 1), screen = 1)
split.screen(c(length(data_bedgraph_variable), 1))
screen(1)

par(mar=c(0.1,1,2,0.1))

i=1
j=1

for (bedg_i in seq_along(data_bedgraph_variable)) {
    screen( i )    
    par(mar=c(0.1, 2, 0.1, 2))
    plotBed(data_bedgraph_variable[[bedg_i]], chrom, chromstart, chromend, 
            row='supplied',  
            color= data_bedgraph_variable[[bedg_i]]$color)
    
    mtext(j, side=2, cex=size_lab)
    i=i+1
    j=j+1
}

dev.off()

