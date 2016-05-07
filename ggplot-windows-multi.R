#********************************************************
#   ggplot-windows-multi.R -- Plot number of each windows along chromosomes using ggplot, support multi-sample input.
#                          
#
#   Author: Nowind
#   Created: 2013-07-10
#   Updated: 2014-11-29
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 14/11/29: The initial version.
#*********************************************************


library(lsr)
library(reshape2)
library(ggplot2)
library(grid)



#*********************************************************
## plot distributions of markers
varfile     <- commandArgs(TRUE)[1]
out_tiff    <- commandArgs(TRUE)[2]

variants <- as.data.frame(read.table(varfile, header=T, sep="\t"))

chrom_num <- length(unique(variants$chrom))
out_height <- chrom_num * 600
tiff(out_tiff, width=4800, height=out_height, units='px', res=600, compression='lzw', bg='transparent')

variants_melted <- melt(variants, id.vars=c("chrom", "interval", "bin"))

ggplot(variants_melted, aes(x=bin, y=value, colour=variable, group=variable)) +
    geom_line() +
    xlab("Windows(x500kbp)") +
    ylab("Number of Variants") +
    scale_colour_hue(name="variable",
                     breaks=unique(variants_melted$variable),
                     labels=unique(variants_melted$variable),
                     l=40) +
    theme_bw() +
    theme(legend.justification=c(1,0), legend.position=c(1,0)) +
    facet_grid(chrom ~ ., scales="free_y")
#*********************************************************

warnings()

