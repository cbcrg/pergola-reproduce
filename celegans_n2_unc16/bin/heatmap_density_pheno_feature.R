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
#############################################################################
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. May 2016                   ###
#############################################################################
### Mean values for each group of worms                                   ###
### Using bed files raw data intercepted with the motion                  ### 
#############################################################################

##Getting HOME directory 
home <- Sys.getenv("HOME")

### Execution example
## Rscript density_heatmaps_tracking.R --files_str1="list_files_str1" --files_str2="list_files_str2"
library(ggplot2)
library("Biostrings")
library("devtools")
library("Gviz")

## bedfiles to GRanges
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
      --path_str1=path_files_str1      - character  
      --path_str2=path_files_str2      - character
      --file_f_str1=file_forward_str1  - character
      --file_f_str2=file_forward_str2  - character
      --file_b_str1=file_backward_str1 - character
      --file_b_str2=file_bacward_str2  - character
      --file_p_str1=file_paused_str1   - character
      --file_p_str2=file_paused_str2   - character
      --help                           - print this text
      
      Example:
      ./plot_speed_motion_mean.R --path_str1=\"path_to_files_str1\" --path_str2=\"list_to_files_str2\" \n")
  
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
  if (is.null (argsL$path_str1)) 
  {
    stop ("[FATAL]: path_str1 arg is mandatory")
  }
  else
  {
    path_str1 <- argsL$path_str1
  }
}

{
  if (is.null (argsL$path_str2)) 
  {
    stop ("[FATAL]: path_str2 not provided", stderr())
  }
  else
  {
    path_str2 <- argsL$path_str2
  }
}

{
  if (is.null (argsL$file_f_str1)) 
  {
    stop ("[FATAL]: file_f_str1 arg is mandatory")
  }
  else
  {
    file_f_str1 <- argsL$file_f_str1
  }
}

{
  if (is.null (argsL$file_b_str2)) 
  {
    stop ("[FATAL]: file_f_str1 arg is mandatory")
  }
  else
  {
    file_f_str2 <- argsL$file_f_str2
  }
}

{
  if (is.null (argsL$file_b_str1)) 
  {
    stop ("[FATAL]: file_b_str1 arg is mandatory")
  }
  else
  {
    file_b_str1 <- argsL$file_b_str1
  }
}

{
  if (is.null (argsL$file_b_str2)) 
  {
    stop ("[FATAL]: file_b_str1 arg is mandatory")
  }
  else
  {
    file_b_str2 <- argsL$file_b_str2
  }
}

{
  if (is.null (argsL$file_p_str1)) 
  {
    stop ("[FATAL]: file_p_str1 arg is mandatory")
  }
  else
  {
    file_p_str1 <- argsL$file_p_str1
  }
}

{
  if (is.null (argsL$file_p_str2)) 
  {
    stop ("[FATAL]: file_p_str1 arg is mandatory")
  }
  else
  {
    file_p_str2 <- argsL$file_p_str2
  }
}

label_gr_2 <- "UNC-16"
label_gr_2 <- "N2"
 
#############################
## Read files bedGraph files
ctrl_worms.base <- path_str2
ctrl.worms.bedg.files <- list.files(ctrl_worms.base, pattern="*.bedGraph$", full.names=TRUE)

scores_bedgr_n2 <- sapply (ctrl.worms.bedg.files, y <- function (x) { data <- read.table (x)                                                      
                                                                      return (data$V4)
})

## They do not have the same number of intervals that is why it doesn't work!!!!
# files do need to have the same number of windows
df_bedg_n2 <- data.frame(scores_bedgr_n2)
df_bedg_n2 [df_bedg_n2 == -10000] <- 0
names(df_bedg_n2) <- paste ("n2_", seq (1 : length(df_bedg_n2[1,])), sep="")
head (df_bedg_n2)
n2_bedg_GR <- GRanges()
n2_bedg_GR <- import(ctrl.worms.bedg.files [[1]], format = "bedGraph")
mcols(n2_bedg_GR) <- df_bedg_n2

## Heatmap
n2_bedg_dt_heatMap <- DataTrack(n2_bedg_GR, name = "midbody speed (microns/s)", type="heatmap", ylim=c(-400, 400), 
                                gradient=c('blue', 'white','red'),
                                #groups = group_tr_rev, col=c(col_case, col_ctrl), 
                                legend=FALSE)

# Read unc16 data
case_worms.base <- path_str1
case.worms.bedg.files <- list.files(case_worms.base, pattern="*.bedGraph$", full.names=TRUE)

scores_bedgr_unc16 <- sapply (case.worms.bedg.files, y <- function (x) { data <- read.table (x)                                                      
                                                                         return (data$V4)
})

df_bedg_unc16 <- data.frame(scores_bedgr_unc16)
df_bedg_unc16 [df_bedg_unc16 == -10000] <- 0
names(df_bedg_unc16) <- paste ("unc16_", seq (1 : length(df_bedg_unc16[1,])), sep="")
unc16_bedg_GR <- GRanges()
unc16_bedg_GR <- import(case.worms.bedg.files[[1]], format = "bedGraph")
mcols(unc16_bedg_GR) <- df_bedg_unc16

## Heatmap unc16
unc16_bedg_dt_heatMap <- DataTrack(unc16_bedg_GR, name = "midbody speed (microns/s)", type="heatmap", ylim=c(-400, 400), 
                                   gradient=c('blue', 'white','red'),
                                   #groups = group_tr_rev, col=c(col_case, col_ctrl), 
                                   legend=FALSE)

### Plot track
# plotTracks(unc16_bedg_dt_heatMap, from=1, to=20000,
#            showSampleNames = TRUE, cex.sampleNames = 0.6)# sample names in heatmap

## join both n2 and unc-16
n2_unc16_bedg_GR <- GRanges()
n2_unc16_bedg_GR <- import(ctrl.worms.bedg.files [[1]], format = "bedGraph")
df_joined <- cbind(df_bedg_n2, df_bedg_unc16)
df_joined [df_joined == -10000] <- 0

mcols(n2_unc16_bedg_GR) <- df_joined

## Heatmap
n2_unc16_bedg_dt_heatMap <- DataTrack(n2_unc16_bedg_GR, name = "midbody speed (microns/s)", 
                                      type="heatmap", ylim=c(-400, 400),
                                      gradient=c('blue', 'white','red'))

gtrack <- GenomeAxisTrack()

read_bed <- function (bed_file) {
  info = file.info(bed_file)
  if (info$size == 0) {
    df_bed <- data.frame (chr="chr1", start=0, end=0, data_type=0, value=0, strand=0, s=0, e=0, color_code=0, body_part="", motion=0, strain="")
  }
  else { df_bed <- read.csv(file=bed_file, header=F, sep="\t")        
         colnames (df_bed) <- c("chr", "start", "end", "data_type", "value", "strand", "s", "e", "color_code",  "body_part", "motion", "strain")
  }
  
  ## We remove this fake rows they were included just to avoid last line of code above to crash
  df_bed <- df_bed [!(df_bed$start == 0 & df_bed$end == 0), ]
  
  return (df_bed)
}

df_f_str1 <- read_bed (file_f_str1)
df_f_str2 <- read_bed (file_f_str2)
df_b_str1 <- read_bed (file_b_str1)
df_b_str2 <- read_bed (file_b_str2)
df_p_str1 <- read_bed (file_p_str1)
df_p_str2 <- read_bed (file_p_str2)

df_f_bed <- rbind (df_f_str1, df_f_str2)
df_b_bed <- rbind (df_b_str1, df_b_str2)
df_p_bed <- rbind (df_p_str1, df_p_str2)

df_f_bed_filt <- df_f_bed [df_f_bed$value != -10000, ]
df_b_bed_filt <- df_b_bed [df_b_bed$value != -10000, ]
df_p_bed_filt <- df_p_bed [df_p_bed$value != -10000, ]

plot_density <- function (df_str1_str2_filt, dir="for") {
  size_titles <- 11
  size_axis <- 8
  size_axis_ticks <- 8
  size_axis_ticks_y <- 8
  size_legend <- 8
  
  cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  labs_plot <- as.vector(levels(df_str1_str2_filt$strain))
 
  if (dir == "back") {
    xmin <- -700
    xmax <- 30 
    tick_interval <- 400
  }
  else if (dir == "for") {
    xmin <- -100
    xmax <- 600
    tick_interval <- 400
  }
  else if (dir == "paused") {
    xmin <- -200
    xmax <- 200
    tick_interval <- 200
  }
  else {
    xmin <- -100
    xmax <- 600
    tick_interval <- 400
  }
  
  breaks_v <- c(-rev(seq(0,abs(xmin), by=tick_interval)[0:-1]), seq (0, xmax, by=tick_interval))
  
  units_lab <- expression(paste(mu, "m / s", "\n", sep=""))
  
  density_str1_str2 <- ggplot(df_str1_str2_filt, aes(x=value, fill=strain)) + geom_density(alpha=0.25) +
    scale_x_continuous (breaks=breaks_v, limits=c(xmin, xmax)) +
    scale_y_continuous(breaks=NULL) +     
    labs (title = "") +
    labs (x = units_lab, y = expression(paste("Probability (", Sigma, "P(x) = 1)", "\n", sep=""))) +         
    theme (axis.text.x = element_text(size=size_axis_ticks)) +
    theme (plot.title = element_text(size=size_titles)) + 
    theme (axis.title.x = element_text(size=size_axis)) +
    theme (axis.title.y = element_text(size=size_axis)) +
    theme (axis.text.x = element_text(size=size_axis_ticks)) +  
    theme (axis.text.y = element_text(size=size_axis_ticks_y)) +
    theme (legend.text = element_text(size=size_legend)) +
    theme (legend.key.size = unit(0.2, "cm"), legend.title = element_blank()) +
    scale_fill_manual(name='',  values = cbb_palette) #+  #labels = labs_plot,
  #   theme (legend.text = element_blank())

  return (density_str1_str2)
}

plot_str1_str2_for <- plot_density (df_f_bed_filt, dir="for")
plot_str1_str2_paused <- plot_density (df_p_bed_filt, dir="paused")
plot_str1_str2_back <- plot_density (df_b_bed_filt, dir="back")

###### Plotting
tiff("heatmap_str1_str2.tiff", width = 18 , height = 14, units = "cm", res=300)

nrows <- 3
ncols <- 3
grid.newpage() 
pushViewport(viewport(layout=grid.layout(nrows, ncols)))
pushViewport(viewport(layout.pos.col=1:3, 
                      layout.pos.row=1:2))

pt <- plotTracks(c(n2_unc16_bedg_dt_heatMap, gtrack), from=1, to=23000, col = NULL,
                 #                  # yellow background
                 #                  background.title = "#4d4d4d", background.panel = "#FFFEDB",
                 showSampleNames = TRUE, cex.sampleNames = 0.6, add=TRUE) ## ADD TRUE!!!! #, title.width=10

popViewport(1)
print(plot_str1_str2_for, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))

print(plot_str1_str2_paused, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

print(plot_str1_str2_back, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))

dev.off()
