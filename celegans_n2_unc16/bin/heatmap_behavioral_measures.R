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
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. July 2017                  ###
#############################################################################
### For each mat files we obtained the correlations between 2 phenotypic. ###
### We obtain the mean along the difference mat files and make a heatmap  ### 
### of the correlation between the difference body parts.                 ###
#############################################################################

#Reading arguments
args <- commandArgs (TRUE) #if not it doesn't start to count correctly

## Default setting when no arguments passed
if ( length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      heatmap_behavioral_corr.R
      
      Arguments:
      --measures_tbl=corr_tbl.csv  - character  
      --help                       - print this text
      
      Example:
      ./heatmap_behavioral_measures.R --corr_tbl=\"path_to_file\" \n")
  
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
  if (is.null (argsL$corr_tbl)) 
  {
    stop ("[FATAL]: corr_tbl arg is mandatory")
  }
  else
  {
    corr_tbl <- argsL$corr_tbl
  }
}

library(ggplot2)
# install.packages("reshape")
# library(reshape)

# Loading params plot:
source("https://raw.githubusercontent.com/cbcrg/mwm/master/lib/R/plot_param_public.R")

corr_data <- read.csv(corr_tbl, sep="\t", header=FALSE, stringsAsFactors = F)

#dcast(corr_data,  V2 ~ V3, fun=mean,drop=FALSE, fill=0)
mean_corr <- aggregate(V4 ~ V2 + V3, data = corr_data, FUN= "mean")
unique_var <- unique(c(mean_corr$V2, mean_corr$V3))
# diagonal <-as.data.frame(cbind(unique_var, unique_var, 1))
#colnames(diagonal) <- colnames(mean_corr)
#mean_corr <- rbind(mean_corr, diagonal)

dpi <- 300
plot_width <- 12
plot_height <- 12 
# name_file <- "heatmap_corr_behavioral_measures"
name_file <- basename(corr_tbl)
file_format <- ".pdf"

name_out <- paste(name_file, file_format, sep="")

ggplot(data = mean_corr, aes(V3, V2, fill = V4, label=V4), color=black, size=12) +
  geom_tile(color = "white")+
  geom_text(aes(label = round(V4, 3))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  labs(x = "", y = "") +
  #theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 20, hjust = 1),
        axis.text.y = element_text(size = 20)) +
  coord_fixed()

ggsave (file=name_out, width = plot_width, height=plot_height)

#corr.m <- dcast(mean_corr, V2 ~ V3)
#row.names(corr.m) <- corr.m [,1]
#corr.m <- corr.m [,-1]
#
#corr.m[!upper.tri(corr.m)] <- 1
#corr.m <- corr.m + t(corr.m) -1
#
#reorder_cor_mat <- function(corr_mat){
#  # Use correlation between variables as distance
#  dd <- as.dist((1-corr_mat)/2)
#  hc <- hclust(dd)
#  corr_mat <-corr_mat[hc$order, hc$order]
#  return (corr_mat)
#}
#
#corr.m <- reorder_cor_mat (corr.m)
#
#upper_triang <- corr.m 
#upper_triang[lower.tri(upper_triang)] <- NA
#upper_tri
## melted_corr.m <- melt(as.matrix(upper_triang), na.rm = TRUE)
#melted_corr.m <- melt(as.matrix(upper_triang))
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
#
#melted_corr.m$value <- as.numeric(as.character(melted_corr.m$value))
#
#melted_corr.m  <- subset (melted_corr.m, value != "NA") 
## Create a ggheatmap
#ggheatmap <- ggplot(melted_cormat, aes(X1, X2, fill = value))+
#  geom_tile(color = "white")+
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
#                       name="Pearson\nCorrelation") +
#  theme_minimal()+ # minimal theme
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                   size = 12, hjust = 1))+
#  coord_fixed()
## Print the heatmap
#print(ggheatmap)
#
################################################
#cormat <- round(cor(mydata),2)
#reorder_cormat <- function(cormat){
#  # Use correlation between variables as distance
#  dd <- as.dist((1-cormat)/2)
#  hc <- hclust(dd)
#  cormat <-cormat[hc$order, hc$order]
#}
#cormat <- reorder_cormat(cormat)
#upper_tri <- get_upper_tri(cormat)
## Melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)
## Create a ggheatmap
#ggheatmap <- ggplot(melted_cormat, aes(X2, X1, fill = value))+
#  geom_tile(color = "white")+
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                       midpoint = 0, limit = c(-1,1), space = "Lab", 
#                       name="Pearson\nCorrelation") +
#  theme_minimal()+ # minimal theme
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                   size = 12, hjust = 1))+
#  coord_fixed()
## Print the heatmap
#print(ggheatmap)
#
#
#
#library(reshape2)
#data.matrix (corr.df)
#
#mydata <- mtcars[, c(1,3,4,5,6,7)]
#head(mydata)
#cor_mat <- round(cor(mydata),2)
#head(cormat)
#
#dd <- as.dist((1-cormat)/2)
#dd
#reorder_cormat <- function(cor_mat){
#  # Use correlation between variables as distance
#  dd <- as.dist((1-cor_mat)/2)
#  hc <- hclust(dd)
#  cor_mat <-cor_mat[hc$order, hc$order]
#}
#cormat <- reorder_cormat(cormat)
#cormat
#
#get_upper_tri <- function(cormat){
#  cormat[lower.tri(cormat)]<- NA
#  return(cormat)
#}
#upper_tri <- get_upper_tri(cormat)
#class(upper_tri)
#melted_cormat <- melt(upper_tri, na.rm = TRUE)






# pair <- combn (unique(c(as.vector(mean_corr$V3), as.vector(mean_corr$V2))), 2, simplify=T)
# pair <- outer (unique(c(as.vector(mean_corr$V3), as.vector(mean_corr$V2))), "paste",)
# 
# # pair <- expand.grid(unique(c(as.vector(mean_corr$V3), as.vector(mean_corr$V2))), unique(c(as.vector(mean_corr$V3), as.vector(mean_corr$V2))))
# #pair <- combn(unique(c(mean_corr$V3, mean_corr$V2)), 2, simplify = F)
# 
# var_cor <- function(x, body_part = mean_corr) { print (x)}
# var_cor <- function(x, body_part = mean_corr) {
#   y <- subset (body_part, V2==x[1] & V3==x[2]) 
#   
#   if (length(y$V2) ==  0) {
#     y <- subset (body_part, V2==x[2] & V3==x[1])
#   }
#   
#   return (y)
# }
# # subset (mean_corr, V2=="tailTip" & V3=="tailTip") 
# # var_cor <- function(x){print(x[3])}
# corr.df <- do.call ( rbind, lapply(pair, var_cor))
# corr.df$V4 <- as.numeric(corr.df$V4)
# 
# # melt(corr.df, na.rm = TRUE)
# 
# dcast(corr.df, V2 ~ V3)
