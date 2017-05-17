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

## Loading libraries
library("ggplot2")
library("ggrepel")

fc_pvalue <- read.table(path2file, header=FALSE)

colnames (fc_pvalue) <- c("variable", "log2FoldChange", "pvalue")
max_y <- max (fc_pvalue$log2FoldChange)

fc_pvalue$highlight <- ifelse(abs(fc_pvalue$log2FoldChange) > 0.2 & -log10(fc_pvalue$pvalue) > 40, "black", "red") 

## Filtering extreme values for plotting name
interesting_p_fc_pos <- subset (fc_pvalue, log2FoldChange > 0.3)
interesting_p_fc_neg <- subset (fc_pvalue, log2FoldChange < -0.18)
interesting_p_pv <- subset (fc_pvalue, -log10(pvalue) > 50)

# png(paste("volcano_plot", ".png", sep=""))
pdf ( paste("volcano_plot", ".pdf", sep=""), height=10, width=12)

with(fc_pvalue, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# with(interesting_p_fc_pos, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
# with(interesting_p_fc_neg, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
# with(interesting_p_pv, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
dev.off()

# png(paste("volcano_plot_pos", ".png", sep=""))
pdf ( paste("volcano_plot_pos", ".pdf", sep=""), height=10, width=12)

if (nrow(interesting_p_fc_pos) != 0) {
with(interesting_p_fc_pos, plot(log2FoldChange, -log10(pvalue), xlim=c(0.4, max_y), pch=20))
with(interesting_p_fc_pos, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
dev.off()
}

volcano_ggplot_pos <- ggplot(interesting_p_fc_pos) +                  
                      geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=highlight)) +                  
                      scale_color_manual(values=c('red','black')) +
                      geom_text_repel(aes(log2FoldChange, -log10(pvalue), label = variable)) +
                      theme_classic(base_size = 16) +
                      theme(legend.position="none")
ggsave (volcano_ggplot_pos, file=paste("volcano_plot_labels_pos", ".png", sep=""))

# png(paste("volcano_plot_neg", ".png", sep=""))
pdf(paste("volcano_plot_neg", ".pdf", sep=""), height=10, width=12)

if (nrow(interesting_p_fc_neg) != 0) {
with(interesting_p_fc_neg, plot(log2FoldChange, -log10(pvalue), xlim=c(-0.42,-0.18), pch=20))
with(interesting_p_fc_neg, text(log2FoldChange, -log10(pvalue), variable, cex=0.8, pos=4, col="red"))
dev.off()
}

volcano_ggplot_neg <- ggplot(interesting_p_fc_neg) +                  
  geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=highlight)) +                  
  scale_color_manual(values=c('red','black'))+
  geom_text_repel(aes(log2FoldChange, -log10(pvalue), label = variable)) +
  theme_classic(base_size = 16) +
  theme(legend.position="none")
# ggsave (volcano_ggplot_neg, file=paste("volcano_plot_labels_neg", ".png", sep=""))
ggsave (volcano_ggplot_neg, file=paste("volcano_plot_labels_neg", ".pdf", sep=""))

velmag <- fc_pvalue [fc_pvalue$variable == "velmag",]

# fc_pvalue_int <- fc_pvalue [abs(fc_pvalue$log2FoldChange) > 0.2 & -log10(fc_pvalue$pvalue) > 40, ]
# fc_pvalue_int <- fc_pvalue [abs(fc_pvalue$log2FoldChange) > 0.8 & -log10(fc_pvalue$pvalue) > 40, ]
fc_pvalue_int <- fc_pvalue [abs(fc_pvalue$log2FoldChange) > 0.4 & -log10(fc_pvalue$pvalue) > 55, ]

volcano_ggplot <- ggplot(fc_pvalue) +                  
                  geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=highlight)) +                  
                  scale_color_manual(values=c('red','black'))+
#                   geom_text_repel(aes(log2FoldChange, -log10(pvalue), label = variable)) +
                  geom_text_repel(data=fc_pvalue_int, aes(log2FoldChange, -log10(pvalue), label = variable)) +
#                   geom_text_repel(data=interesting_p_fc_pos, aes(log2FoldChange, -log10(pvalue), label = variable)) +
#                   geom_text_repel(data=velmag, aes(log2FoldChange, -log10(pvalue), label = variable)) +
                  theme_classic(base_size = 16) +
                  theme(legend.position="none") +
#                   xlim(c(-1,1.6))
                  xlim(c(-0.5,1.5))

# ggsave (volcano_ggplot, file=paste("volcano_plot_labels", ".png", sep=""))
ggsave (volcano_ggplot, file=paste("volcano_plot_labels", ".pdf", sep=""))

fc_pvalue_int <- fc_pvalue [abs(fc_pvalue$log2FoldChange) > 0.8 & -log10(fc_pvalue$pvalue) > 80, ]

volcano_ggplot <- ggplot(fc_pvalue) +                  
  geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=highlight)) +                  
  scale_color_manual(values=c('red','black'))+
  geom_text_repel(data=fc_pvalue_int, aes(log2FoldChange, -log10(pvalue), label = variable)) +
  theme_classic(base_size = 16) +
  theme(legend.position="none") +
  xlim(c(-0.5,1.5))

ggsave (volcano_ggplot, file=paste("volcano_plot_stringent_threshold", ".pdf", sep=""))

fc_pvalue$pvalue <- -log10(fc_pvalue$pvalue)
write.table(fc_pvalue, "tbl_fc_pvalues.txt", sep="\t", col.names=TRUE, row.names=FALSE) 

