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
### t-test between scores falling in annotated jaaba regions vs those    ###
### not annotated                                                        ### 
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
        ttest_var_annotated_jaaba
        
        Arguments:
        --path2files=someValue     - character, path to read bedGraph files        
        --variable_name=someValue  - character, variable name
        --help                     - print this text
        
        Example:
        ./ttest_var_annotated_jaaba.R --path2files=\"/foo/variables\" --variable_name=\"var_name\"\n")
    
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
    if (is.null (argsL$path2files)) 
    {
        stop ("[FATAL]: Path to files variables bedgraph files is mandatory")
    }
    else
    {
        path2files <- argsL$path2files
    }
}

# variable name
{
    if (is.null (argsL$variable_name)) 
    {
        print ("[WARNING]: Variable to plot is set to \"default=Speed\"\n")
        variable_name <- "Speed"
    }
    else
    {        
        variable_name <- argsL$variable_name
    }
}

## Loading libraries
library("ggplot2")

## loading aes parameters for plotting 
source("https://raw.githubusercontent.com/cbcrg/mwm/master/lib/R/plot_param_public.R")

# setwd(path2files)

## this are the two in principle that should be different
# variable <- "velmag"
# variable <- "dtheta"
variable <- variable_name

## I should take into account also the tracks without any annotation
## as for instances bed 1

# variable <- "velmag_tail"
# variable <- "dtheta"
# variable <- "phi"
# variable <- "dell2nose"
# variable <- "angle2wall"
# path2files <- "/Users/jespinosa/2017_sushi_pergola/lib/python"
# path2files <- "/Users/jespinosa/2017_sushi_pergola/lib/nxf/work/e0/2887dcdce79d24a889eac82aa20699/results_annot"
files_annotated <- list.files(path=path2files, pattern=paste("values_.*", variable, ".txt$", sep=""), full.names = TRUE)

data_bed <- lapply(files_annotated, read.csv, header=FALSE, sep="\t", stringsAsFactors=FALSE)
data_bed.df <- do.call(rbind, data_bed)

v_annotated <- abs(as.numeric(unlist(strsplit(data_bed.df$V10, ","))))

files_annotated_comp <- list.files(path=path2files, pattern=paste("values_.*", variable, ".comp.txt$", sep=""), full.names = TRUE)
data_bed_comp <- lapply(files_annotated_comp, read.csv, header=FALSE, sep="\t", stringsAsFactors=FALSE)

data_bed_comp.df <- do.call(rbind, data_bed_comp)

v_no_annotated <- abs(as.numeric(unlist(strsplit(data_bed_comp.df$V4, ","))))

t_result <- t.test(v_annotated, v_no_annotated)

## Fold change calculation
# to avoid zeros I substitute 0 values by the mean of the values
v_annotated[v_annotated==0] <- mean(v_annotated)
v_no_annotated[v_no_annotated==0] <- mean(v_no_annotated)

log2FoldChange = mean(log2(v_annotated)) - mean(log2(v_no_annotated))
vector_FC <- c(variable, log2FoldChange, t_result$p.value)
M_t_vector_FC <- as.matrix(t(vector_FC))
write.table(M_t_vector_FC, stdout(), sep="\t", col.names=FALSE, row.names=FALSE)

## bar plot group comparison
group <- c(rep("Annotated", length (v_annotated)), rep("No annotated", length (v_no_annotated)))
df_values <- data.frame(id = group, value = c(v_annotated, v_no_annotated))

## outliers out 
ylim1 = boxplot.stats(df_values$value)$stats[c(1, 5)]

## colors
cbb_palette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
name_out <- paste (variable, ".", "png", sep="")

ggplot(df_values, aes(id, value, fill=id)) + geom_boxplot() + 
    labs (title = "Jaaba annotated vs. non-annotated intervals\n", 
          y = paste(variable, "\n", sep=""), x = "\nGroup") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = cbb_palette) +
    theme(legend.position="none") +
    # outliers out
    coord_cartesian(ylim = ylim1*1.05) +
    annotate("text", x=2.3, y=ylim1[2], label=paste("p-value=", signif (t_result$p.value,3)))

ggsave (file=name_out)