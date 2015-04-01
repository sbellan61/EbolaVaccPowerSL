####################################################################################################
## Plot Figure S7.
####################################################################################################
## Code base accompanying:
## 
## Bellan, SE, JRC Pulliam, CAB Pearson, DChampredon, SJ Fox, L Skrip, AP Galvani, M Gambhir, BA
## Lopman, TC Porco, LA Meyers, J Dushoff (2015). The statistical power and validity of Ebola
## vaccine trials in Sierra Leone: A simulation study of trial design and analysis. _Lancet
## Infectious Diseases_.
##
## Steve Bellan, March 2015
## License at bottom.
####################################################################################################
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)

labs <- c('','log')
thing <- 'FalsePosFluct'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))
source('ggplotTheme.R')

pf$cvWeeklyDecay <- pf[, as.numeric(levels(cvWeeklyDecay)[cvWeeklyDecay])]
pf$weeklyDecay <- pf[, as.numeric(levels(weeklyDecay)[weeklyDecay])]
pf$cvClus <- pf[, as.numeric(levels(cvClus)[cvClus])]
pf$cvClusTime <- pf[, as.numeric(levels(cvClusTime)[cvClusTime])]

####################################################################################################
## Figure S7 - Type I errors by phenomenological amount of fluctuation False positive rates for
## varying degrees of inter- and intra-cluster heterogeneity in trends.
####################################################################################################
subs <- pf[,trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], 
                aes(cvWeeklyDecay, stoppedNAR, colour=trial, linetype=order)) + thsb +
    scale_x_continuous(limits=c(0,1), minor_breaks=NULL, breaks = unique(pf[,cvWeeklyDecay])) +  
    xlab('coef variation of weekly decay rate') + ylab('False Positive Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(mod ~ cvClusTime, scales='free') + scale_color_manual(values=group.colors)
 p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
ggsave(paste0('Figures/Fig S7 - Type I by fluct & weekly decay CV.pdf'), p.tmp, w = 7.5, h = 5.5)
ggsave(paste0('Figures/Fig S7 - Type I by fluct & weekly decay CV.png'), p.tmp, w = 7.5, h = 5.5)


####################################################################################################
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/sbellan61/EbolaVaccPowerSL
####################################################################################################
