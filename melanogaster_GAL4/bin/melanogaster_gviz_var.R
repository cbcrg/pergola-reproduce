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
      --help                               - print this text
      
      Example:
      ./melanogaster_gviz_var.R --path_bed_files=\"path_bed_files\" \n")      
      
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
  tr <- DataTrack(bedgraph_GR, type="a",
            #background.title = l_gr_color[[i_group_exp]],           
            showAxis = F, name = id)
  return (tr)
} )

names(bedgraph_tracks) <- as.numeric(gsub(".+tr_(\\d+)(_.+$)", "\\1", bedgraph.var_files))
bedgraph_tracks <- bedgraph_tracks[as.character(1:length(bedgraph.var_files))]

name_file <- "gviz_jaaba_var.tiff"
tiff(name_file, width = 45 , height = 34, units = "cm", res=300)
plotTracks(bedgraph_tracks, stacking="dense", collapse=FALSE, shape = "box", col="red")
dev.off()

#######################
### Lo que habria que hacer es en la fucnion leerlas a la vez y entonces dentro de ella crear la data track 
### con las dos group lo que sea

# path2bedg_files <- "/Users/jespinosa/scratch/ba/c0bb5e7aac265be4a11d77cba09e9f/results_bedGr"
# # path2bedg_files <- "/Users/jespinosa/scratch/94/09f279eec46db44d013433941f7326/results_var"
# variable_name <- "velmag"
# base_folder <- path2bedg_files
# # chase.bedgraph.variable.files <- mixedsort(list.files(base_folder, pattern=paste("values.*", variable_name, ".bedGraph$", sep=""), full.names=TRUE))
# chase.bedgraph.variable.files.comp <- mixedsort(list.files(base_folder, pattern=paste("values.*", variable_name, "*.comp.*.bedGraph$", sep=""), full.names=TRUE))
# 
# ### TEST comp tracks, those no annotated as chasing
# ## the rest are the ones annotated as chasing
# # chase.bedgraph.variable.files <- mixedsort(list.files(base_folder, pattern=paste("values.*", variable_name, ".bedGraph$", sep=""), full.names=TRUE))
# 
# # bedgraph_tracks_comp <- lapply(chase.bedgraph.variable.files.comp, function (bedgraph) {
# #   #   id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
# #   #   name_id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedg)))
# #   id <- gsub("values_", "", gsub(paste("_", variable_name, ".comp.bedGraph", sep=""), "", basename(bedgraph)))
# #   #   id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedgraph)))
# #   print (id)
# #   bedgraph_GR <- import(bedgraph, format = "bedGraph")
# #   tr <- DataTrack(bedgraph_GR, type="a", #ylim = c(min_heatmap, max_heatmap),
# #                   #background.title = l_gr_color[[i_group_exp]],
# #                   #gradient=c(color_min, color_max), 
# #                   showAxis = F, name = id)
# #   
# #   #   tr <- AnnotationTrack(bedgraph_GR, name = paste ("", id, sep="")) #,
# #   # #                         fill=bed_GR$itemRgb)
# #   #                       background.title = c)  
# #   #   names(tr) <- id
# #   return (tr)
# # } )
# 
# # names(bedgraph_tracks_comp) <- as.numeric(gsub("values_", "", gsub(paste("_", variable_name, ".comp.bedGraph", sep=""), "", basename(chase.bedgraph.variable.files.comp))))
# # bedgraph_tracks_comp <- bedgraph_tracks_comp[as.character(1:length(bedgraph.var_files))]
# # 
# # # collapse may allow to over the tracks with the same id 
# # plotTracks(c(bedgraph_tracks, bedgraph_tracks_comp), stacking="dense", collapse=TRUE, shape = "box", col="red")
# # plotTracks(c(bedgraph_tracks, bedgraph_tracks_comp), from=0, to=1000)
# # plotTracks(c(bedgraph_tracks_comp), from=0, to=1000)
# # plotTracks(c(bedgraph_tracks), from=0, to=1000)
# 
# bedgraph_tracks_annot <- lapply(chase.bedgraph.variable.files, function (bedgraph) {
#     #   id <- gsub(".+tr_(\\d+)(_.+$)", "\\1", bed)
#     #   name_id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedg)))
#     id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedgraph)))
#     #   id <- gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(bedgraph)))
#     print (id)
#     bedgraph_GR <- import(bedgraph, format = "bedGraph")
#     tr <- DataTrack(bedgraph_GR, type="a", #ylim = c(min_heatmap, max_heatmap),
#                     #background.title = l_gr_color[[i_group_exp]],
#                     #gradient=c(color_min, color_max), 
#                     showAxis = F, name = id)
#     
#     #   tr <- AnnotationTrack(bedgraph_GR, name = paste ("", id, sep="")) #,
#     # #                         fill=bed_GR$itemRgb)
#     #                       background.title = c)  
#     #   names(tr) <- id
#     return (tr)
# } )
# 
# names(bedgraph_tracks_annot) <- as.numeric(gsub("values_", "", gsub(paste("_", variable_name, ".bedGraph", sep=""), "", basename(chase.bedgraph.variable.files))))
# bedgraph_tracks_annot <- bedgraph_tracks_annot[as.character(1:length(bedgraph.var_files))]
# 
# # collapse may allow to over the tracks with the same id 
# plotTracks(c(bedgraph_tracks, bedgraph_tracks_annot), stacking="dense", collapse=TRUE, shape = "box", col="red")
# plotTracks(c(bedgraph_tracks, bedgraph_tracks_annot), from=0, to=1000)
# plotTracks(c(bedgraph_tracks_annot), from=0, to=5000)
# 
# plotTracks(c(bedgraph_tracks, bedgraph_tracks_annot), 
#            #            groups = rep(c("control", "treated"), each = 20), 
#            from=0, to=5000)
# plotTracks(c(bedgraph_tracks), from=0, to=1000)
# 
