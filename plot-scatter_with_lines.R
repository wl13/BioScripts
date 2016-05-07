#********************************************************
#   ggplot-windows-count.R -- Plot number of varaints in each window along each chromosomes using ggplot.
#                          
#
#   Author: Nowind
#   Created: 2013-07-10
#   Updated: 2016-01-03
#   Version: 1.1.2
#
#   Change logs:
#   Version 1.0.0 14/11/28: The initial version.
#   Version 1.1.0 15/12/30: Updated: change input/output file handle.
#   Version 1.1.1 16/01/02: Updated: baseline is now read in from commandline.
#   Version 1.1.2 16/01/03: Updated: change some colors.
#*********************************************************


library(lsr)
library(reshape2)
library(ggplot2)
library(grid)


#*********************************************************
## plot distributions of markers
infile   <- file("stdin")
baseline <- as.numeric(commandArgs(TRUE)[1])
outfile  <- commandArgs(TRUE)[2]


dat <- as.data.frame(read.table(infile, header=T, sep="\t"))



ggplot(dat, aes(x=bin, y=val_mean)) +
    geom_errorbar(data=subset(dat, group == "group1"), aes(ymin=val_mean-val_sem, ymax=val_mean+val_sem), color= "#E41A1C", width=.3, size=0.5) +
    geom_line(data=subset(dat, group == "group1"), color= "#E41A1C", size = 0.5) +
    geom_point(data=subset(dat, group == "group1"), color= "#E41A1C", size=2, shape=21, fill="white") +
    #geom_errorbar(data=subset(dat, group == "group2"), aes(ymin=val_mean-val_sem, ymax=val_mean+val_sem), color= "#07F900", width=.1, size=0.5) +
    #geom_line(data=subset(dat, group == "group2"), color= "#07F900", size = 0.5) +
    #geom_point(data=subset(dat, group == "group2"), color= "#07F900", size=3, shape=21, fill="white") +
    xlab("xlab") + ylab("ylab") +
    scale_y_continuous(breaks=seq(0.01, 0.02, 0.001)) +
    scale_x_continuous(breaks=seq(0, 20, 2)) +
    geom_hline(yintercept=baseline, linetype=5, size=0.5, color="#689EC9") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "grey")) +
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8)) +
    theme(legend.justification=c(1,1),
          legend.position=c(1,1))


ggsave(file=outfile, plot=last_plot(), width=120, height=100, units="mm", dpi=600, bg='transparent')

#*********************************************************


###print(summary(dat))

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right



warnings()

