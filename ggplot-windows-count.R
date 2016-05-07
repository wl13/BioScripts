#********************************************************
#   ggplot-windows-count.R -- Plot number of varaints in each window along each chromosomes using ggplot.
#                          
#
#   Author: Nowind
#   Created: 2013-07-10
#   Updated: 2014-11-28
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 14/11/28: The initial version.
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

chrom_num <- length(unique(variants$CHROM))
out_height <- chrom_num * 600
tiff(out_tiff, width=4800, height=out_height, units='px', res=600, compression='lzw', bg='transparent')

ggplot(variants, aes(x=BIN_ID, y=COUNT, colour=Sample, group=Sample)) +
    geom_line() +
    xlab("Windows(x500kbp)") +
    ylab("Number of Markers") +
    scale_colour_hue(name="Sample",
                     breaks=unique(variants$Sample),
                     labels=unique(variants$Sample),
                     l=40) +
    theme_bw() +
    theme(legend.justification=c(1,0), legend.position=c(1,0)) +
    facet_grid(CHROM ~ .)
#*********************************************************

warnings()

