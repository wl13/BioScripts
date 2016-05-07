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
infile      <- commandArgs(TRUE)[1]
out_tiff    <- commandArgs(TRUE)[2]

dat <- as.data.frame(read.table(infile, header=T, sep="\t"))

out_height <- 3200
tiff(out_tiff, width=3200, height=out_height, units='px', res=600, compression='lzw', bg='transparent')


ggplotRegression <- function (fit) {

    require(ggplot2)

    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point(shape=1) +
        xlab("Divergence") +
        ylab("Diversity") +
        stat_smooth(method = "lm", se=FALSE, col = "red") +
        labs(title = paste(
                 "R2 = ",signif(summary(fit)$r.squared, 3),
                 ", Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
                 ", Intercept =",signif(fit$coef[[1]], 3),
                 ", Slope =",signif(fit$coef[[2]], 3),
                 ", P =",signif(summary(fit)$coef[2,4], 3), sep = "")) +
        theme(plot.title = element_text(size = 10))
}

fit <- lm(Divergence ~ Diversity, data = dat)
ggplotRegression(fit)
#*********************************************************


#lm_eqn = function(df){
#    m = lm(y ~ x, df);
#    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#         list(a = format(coef(m)[1], digits = 2), 
#              b = format(coef(m)[2], digits = 2), 
#             r2 = format(summary(m)$r.squared, digits = 3)))
#    as.character(as.expression(eq));                 
#}
#
#
#
#ggplot(dat, aes(x=x, y=y)) +
#    geom_point(shape=1) +
#    scale_colour_hue(l=50) +
#    xlab("Divergence(1Mbp,all)") +
#    ylab("Diversity(1Mbp,all)") +
#    geom_smooth(method='lm',se=FALSE, color="black", formula=y ~ x) +
#    theme_bw() +
#    theme(legend.justification=c(1,0), legend.position=c(1,0)) +
#    geom_text(aes(x = 0.03, y = 0.0045, label = lm_eqn(dat)), parse = TRUE)

warnings()

