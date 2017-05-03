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
### Volcano plot from comparison between annotated and not annotated     ###
### regions                                                              ### 
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
      volcano_plot_jaaba
      
      Arguments:
      --path2file=someValue     - character, path to read bedGraph files              
      --help                     - print this text
      
      Example:
      ./volcano_plot_jaaba.R --path2file=\"/foo/variables\"\n")
  
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
  if (is.null (argsL$path2file)) 
  {
    stop ("[FATAL]: Path to file containing P-values and fold-changes is mandatory")
  }
  else
  {
    path2file <- argsL$path2file
  }
}

fc_pvalue <- read.table(path2file, header=FALSE)
colnames(fc_pvalue) <- c("variable", "log2FoldChange", "pvalue")

## Filtering extreme values for plotting name
interesting_p_fc_pos <- subset (fc_pvalue, log2FoldChange > 0.3)
interesting_p_fc_neg <- subset (fc_pvalue, log2FoldChange < -0.18)
interesting_p_pv <- subset (fc_pvalue, -log10(pvalue) > 50)

png(paste("volcano_plot", ".png", sep=""))
with(fc_pvalue, plot(log2FoldChange, -log10(pvalue), pch=20, xlim=c(-0.65,0.65), main="Volcano plot"))
# with(interesting_p_fc_pos, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
# with(interesting_p_fc_neg, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
# with(interesting_p_pv, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
dev.off()

png(paste("volcano_plot_pos", ".png", sep=""))
with(interesting_p_fc_pos, plot(log2FoldChange, -log10(pvalue), xlim=c(0.25,0.65), pch=20))
with(interesting_p_fc_pos, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
dev.off()

png(paste("volcano_plot_neg", ".png", sep=""))
with(interesting_p_fc_neg, plot(log2FoldChange, -log10(pvalue), xlim=c(-0.63,-0.18), pch=20))
with(interesting_p_fc_neg, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
dev.off()

fc_pvalue$pvalue <- -log10(fc_pvalue$pvalue)
write.table(fc_pvalue, "tbl_fc_pvalues.txt", sep="\t", col.names=TRUE, row.names=FALSE) 
