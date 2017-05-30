#********************************************************
#   ggplot-windows-count.R -- Plot number of varaints in each window along each chromosomes using ggplot.
#                          
#
#   Author: Nowind
#   Created: 2013-07-10
#   Updated: 2017-05-30
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 14/11/28: The initial version.
#   Version 1.1.0 17/05/30: Add support for other figure format.
#*********************************************************


library(lsr)
library(reshape2)
library(ggplot2)
library(grid)
library(tools)



#*********************************************************
## plot distributions of markers
varfile    <- commandArgs(TRUE)[1]
out_fig    <- commandArgs(TRUE)[2]

variants <- as.data.frame(read.table(varfile, header=T, sep="\t"))

chrom_num <- length(unique(variants$CHROM))
out_height <- chrom_num * 600

if (file_ext(out_fig) == "png") {
    png(out_fig, width=4800, height=out_height, units='px', res=600, bg='transparent')
} else if (file_ext(out_fig) == "tiff") {
    tiff(out_fig, width=4800, height=out_height, units='px', res=600, compression='lzw', bg='transparent')
} else if (file_ext(out_fig) == "jpg") {
    jpeg(out_fig, width=4800, height=out_height, units='px', res=600, bg='white')
} else if (file_ext(out_fig) == "bmp") {
    bmp(out_fig, width=4800, height=out_height, units='px', res=600, bg='transparent')
} else if (file_ext(out_fig) == "svg") {
    svg(out_fig, width=12, height=out_height/400, pointsize=12, onefile=TRUE, family="sans", bg='transparent')
} else if (file_ext(out_fig) == "pdf") {
    cairo_pdf(out_fig, width=12, height=out_height/400, pointsize=12, onefile=TRUE, family="sans", bg='transparent', fallback_resolution = 600)
} else if (file_ext(out_fig) == "ps") {
    cairo_ps(out_fig, width=12, height=out_height/400, pointsize=12, onefile=TRUE, family="sans", bg='transparent', fallback_resolution = 600)
}

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

