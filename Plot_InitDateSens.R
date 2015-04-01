####################################################################################################
## Plot sensitivity analysis to trial start date.
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
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid); library(scales)
thing <- 'initDateSens'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))
source('ggplotTheme.R')
pf[, length(design), list(trial, remProtDel, remStartFin)]

####################################################################################################
## Figure 6
####################################################################################################
subs <- pf[,  immunoDelay==21 & ((trial == 'SWCT' & remProtDel==T & mod=='relabCoxME') | (trial == 'RCT' & order=='risk-prioritized' & mod =='CoxME'))]
p.tmp <- ggplot(pf[subs]) +
  aes(x=trialStartDate, y=vaccGoodNAR, colour=trial, linetype=order) + 
  thsb + theme(axis.text.x = element_text(angle=90)) +
  scale_x_date(labels = date_format("%b-%d"), breaks = pf[,unique(trialStartDate)], minor_breaks=NULL) +
  scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) +  
  xlab('trial start date') + ylab('power') + 
  geom_rect(aes(xmin=as.Date('2015-02-18'), xmax = as.Date('2015-03-18'), ymin=0, ymax=1), fill = "lightgrey", color=NA) +
  geom_line(size=1) +
  scale_color_manual('', values=group.colors) +
      guides(colour = guide_legend(override.aes = list(linetype=c(2,1)))) +
          theme(legend.justification=c(2,1), legend.position=c(1,1.25)) +
  scale_linetype_discrete(guide=F)     
p.tmp
ggsave(paste0('Figures/Fig 6 - Power by start date SL.pdf'), p.tmp, w = 4, h = 3)

pf[subs, list(vaccGoodNAR), list(trial, pit,trialStartDate)]

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
