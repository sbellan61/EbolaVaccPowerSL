####################################################################################################
## Plot power & false positive rates from main analyses.
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
thing <- 'All'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))
source('ggplotTheme.R')
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.025 & mod=='CoxME']

## look at all simulation types (36 of each from 4 vaccine efficacies X 9 models)
arrange(pf[, length(meanNAR), list(trial, remProtDel, remStartFin,propInTrial, ord, delayUnit, immunoDelay)], propInTrial)
nlevels(pf$mod)
unique(pf$vaccEff)

nrow(pf)*2040/9 ## 399,840 simulations, with 9 models each

## Set up line types for ease
pf$analysis <- 1
pf[trial=='SWCT' & remProtDel==T & remStartFin==F, analysis:=1]
pf[trial=='SWCT' & remProtDel==T & remStartFin==T, analysis:=2]
pf[trial=='SWCT' & remProtDel==F & remStartFin==F, analysis:=3]

pf$analysis <- factor(pf$analysis, levels = 1:3, labels = c('(A) SWCT including \nearly/late person-time \n& excluding protective delay'
                                                     ,'(B) SWCT excluding \nboth protective delay \n& early/late person-time'
                                                     ,'(C) SWCT including protective delay \n& early/late person-time'
                                                     ))
## Relabel models
pf$modLab2 <- pf$mod
levels(pf$modLab2)[levels(pf$modLab2) %in% c('CoxME','bootCoxME','relabCoxME')] <- c('(A) CoxPH', '(B) Bootstrap', '(C) Permutation')
## Only include SWCT remPD_SF pt
pf$mainAn <- pf[, trial!='SWCT' | (trial=='SWCT' & as.numeric(analysis)==1) & propInTrial <= .1]

####################################################################################################
## Figure 4 - False positive rates by design & analysis.
####################################################################################################
subs <- pf[, mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=trial, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.025,.1), minor_breaks=NULL, breaks = c(.025,.05,.075,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('False Positive Rate') + 
    #scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~modLab2, scales = "free_y") + scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
p.tmpunt
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Figure 4 - Type I by analysis',typ), p.tmpunt, w = 8.5, h = 3.5)

subs <- pf[, mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','GLMFclus', 'bootCoxME','bootGLMFclus', 'relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
pf[vaccEff==0 & immunoDelay==21 & subs & pit=='10%', list(trial, stoppedNAR,mod)]

####################################################################################################
## Figure S3 - SWCT false positive rates by person-time included
####################################################################################################
subs <- pf[, trial %in% c('SWCT') & vaccEff==0 & mod %in% c('CoxME','GLMFclus','bootCoxME','relabCoxME') & immunoDelay==21]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = percent, limits=c(0,.3), minor_breaks=NULL, breaks = c(0,.1,.2,.3)) +  
    xlab('% of district-level cases in trial population') + ylab('False Positive Rate') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) + theme(legend.key.height=unit(3,"line")) + 
    geom_line(size=1, colour = group.colors['SWCT']) + facet_wrap(~mod, scales = "free_y", nrow = 1) + scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.2))
p.tmpunt
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Figure S3 - type I by SWCT person-time',typ), p.tmpunt, w = 12, h = 3.5)

####################################################################################################
## Type I errors panels by trial type
thax <- element_text(colour = 'black', size = 8)
subs <- pf[, mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
temp.cols <- c('Cox PH Frailty'='orange', 'permutation test over\nCox PH Frailty'='purple',
               'bootstrap over\nCox PH Frailty' = 'dodger blue')
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=model, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.025,.1), minor_breaks=NULL, breaks = c(.025,.05,.075,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~trial, scales = "free_y") + scale_color_manual(values=temp.cols)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
p.tmpunt
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Type I by trial SL', typ), p.tmpunt, w = 6.5, h = 3.5)

####################################################################################################
## Figure 5 - Power
####################################################################################################
subs <- pf[,  mainAn==T & vaccEff>0 &  immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT','FRCT') & mod =='CoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=trial, linetype=order)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit, scales = "free_y",nrow=1) + 
    ggtitle('expected % of district-level cases in trial population') + 
    scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/Figure 5 - Power SL', typ), p.tmpunt, w = 9, h = 4)

tmpS <- arrange(pf[subs & trial=='SWCT', list(trial, vaccEff, vaccGoodNAR, order,pit)],pit)
tmpR <- arrange(pf[subs & trial=='RCT' & order=='risk-prioritized', list(trial, vaccEff, vaccGoodNAR, order,pit)], pit)
range(tmpR$vaccGoodNAR/tmpS$vaccGoodNAR, na.rm=T) ## range of fold increased in power of RCT vs SWCT
## 3-10X

####################################################################################################
## Figure 10 - Power by protective delay
####################################################################################################
for(ii in 1:2) {
    if(ii==1) pdf('Figures/Figure 10 - Power by seroconversion delay SL.pdf', w = 4, h=5)
    if(ii==2) png('Figures/Figure 10 - Power by seroconversion delay SL.png', w = 4, h=5, res = 300, units='in')
    par(mar = c(5,4,0,0))
    subs <- pf[, propInTrial==.05 & vaccEff == .9 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial=='RCT' & order=='risk-prioritized' & mod =='CoxME'))]
    pf[subs]
    plot(0,0, type = 'n', xlim = c(5, 21), ylim = c(0,1), bty='n', las = 1, axes=F,
         xlab ='delay until protective efficacy', ylab='power to detect effective vaccine')
    axis(1, c(5,21))
    axis(2, seq(0,1,l=5),las=1)
    with(pf[subs & trial=='RCT'], lines(immunoDelay, vaccGoodNAR, col = group.colors['RCT'], lty = 2, lwd = 3))
    with(pf[subs & trial=='SWCT' & as.numeric(analysis)==1], lines(immunoDelay, vaccGoodNAR, col = group.colors['SWCT'], lty = 1, lwd = 3))
    with(pf[subs & trial=='SWCT' & as.numeric(analysis)==2], lines(immunoDelay, vaccGoodNAR, col = group.colors['SWCT'], lty = 2, lwd = 3))
    with(pf[subs & trial=='SWCT' & as.numeric(analysis)==3], lines(immunoDelay, vaccGoodNAR, col = group.colors['SWCT'], lty = 3, lwd = 3))
    graphics.off()
}

subs <- pf[, propInTrial==.05 & vaccEff == .9 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial=='RCT' & order=='risk-prioritized' & mod =='CoxME'))]
arrange(pf[subs & (trial=='RCT' | trial=='SWCT' & as.numeric(analysis)==1), list(trial, immunoDelay, propInTrial, vaccGoodNAR)],trial)

####################################################################################################
##  Power by SWCT design all analyses
subs <- pf[, vaccEff>0 & propInTrial <= .1 & immunoDelay==21 & (trial == 'SWCT' & mod %in% c('CoxME','GLMFclus','bootCoxME','relabCoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_line(size=1, colour=group.colors['SWCT']) + facet_wrap(mod~pit, scales = "free",nrow=4) +
 theme(legend.key.height=unit(3,"line")) + 
    ggtitle('expected % of district-level cases in trial population') + 
    scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/Power by SWCT pt', typ), p.tmpunt, w = 9, h = 8)

####################################################################################################
##  Power by SWCT design permutation
subs <- pf[, vaccEff>0 & propInTrial <= .1 & immunoDelay==21 & (trial == 'SWCT' & mod %in% 'relabCoxME')]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_line(size=1, colour=group.colors['SWCT']) + facet_wrap(~pit, scales = "free", nrow = 1) +
 theme(legend.key.height=unit(3,"line")) + 
    ggtitle('expected % of district-level cases in trial population') + 
    scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.43), breaks=seq(0,.4,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/Power by SWCT pt (permutation)', typ), p.tmpunt, w = 9, h = 3.5)

####################################################################################################
##  Prob rej null by SWCT design
subs <- pf[, immunoDelay==21 & propInTrial <= .1 & (trial == 'SWCT' & mod %in% c('CoxME','GLMFclus','bootCoxME','relabCoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(0,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('proportion rejecting null hypothesis that vaccine does not affect risk') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_line(size=1, colour=group.colors['SWCT']) + facet_wrap(mod~pit, scales = "free",nrow=4) +
    theme(legend.key.height=unit(3,"line")) + 
    ggtitle('expected % of district-level cases in trial population') + 
    geom_hline(yintercept=.025, color='dark gray', size = 1) + 
    scale_color_manual(values=group.colors)
p.tmplog <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25, .5, 1), limits = c(.005,1), minor_breaks=NULL) 
p.tmplog
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/CoxPH Prob rej null by SWCT pt log', typ), p.tmplog, w = 9, h = 8)

####################################################################################################
## Figure S5 - Power by model
####################################################################################################
pf$modelLab <- pf$model
pf$class <- 'parametric'
pf[grepl('boot',mod), class:='bootstrap']
pf[grepl('relab',mod), class:='permutation']
pf$class <- factor(pf$class, levels = c("parametric", "bootstrap", "permutation"))
pf$model <- gsub('boot','',pf$mod)
pf$model <- factor(gsub('relab','',pf$model))
levels(pf$model)[1:2] <- c('CoxPH','Poisson regression')
subs <- pf[, mainAn==T & immunoDelay==21 & (trial == 'SWCT' | (trial=='RCT' & order=='risk-prioritized' & class!='bootstrap'))]
subs <- subs & pf[,model %in% c('CoxPH','Poisson regression')]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=class, linetype=model)) +
    thsb +
    scale_x_continuous(labels = formatC, limits=c(0,1), breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('proportion rejecting null hypothesis that vaccine does not affect risk') + 
    scale_linetype_manual(values=1:3) + ## scale_color_manual(values=group.colors) +
    geom_hline(yintercept=.025, color='dark gray', size = 1) + 
    scale_color_manual(values=c(parametric='orange',bootstrap='dodger blue',permutation='purple')) + 
    geom_line(size=.9)  + facet_wrap(trial ~ pit, scales = "free_y", nrow = 2) #,nrow=1) + 
#p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1))
p.tmplog <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25, .5, 1), limits = c(.005,1), minor_breaks=NULL) 
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Figure S5 - Power by model class',typ), p.tmpunt, w = 8, h = 6)
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Figure S5 - Power by model class log',typ), p.tmplog, w = 8, h = 6)

####################################################################################################
## Figure S8 - Show power by # of cases
####################################################################################################
subs <- pf[,  propInTrial<=.1 & mainAn==T & vaccEff==.9 & immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT') & mod =='CoxME'))]
pf$type <- pf[, paste(order,trial)]
pf[grepl('SWCT',type), type:='SWCT']
group.colors2 <- c("orange", 'purple','blue', "#F8766D")
names(group.colors2) <- pf[subs,unique(type)]
p.tmp <- ggplot(pf[subs], aes(caseTot, vaccGoodNAR, colour=type)) + thsb +
    xlab('# of cases in trial population') + ylab('power to detect effective vaccine') + 
    geom_line(size=1.5)  +
    geom_point(aes(caseTot, vaccGoodNAR, colour=type, shape=pit), size = 5) +
    scale_color_manual(values=group.colors2) + scale_shape_manual(values=c(15:18)) + 
    scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) 
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Fig S8 - Power by cases',typ), p.tmp, w = 8, h = 4)

####################################################################################################
## Figure S9 - Cases by trial type.
####################################################################################################
subs <- pf[, mainAn==T & immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT','FRCT') & mod =='CoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, caseTot, colour=trial, linetype=order)) + thsb + 
    scale_x_continuous(labels = formatC, limits=c(0,.9), breaks = pf[,unique(vaccEff)],minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('# EVD cases') +
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit, scale='free_y', nrow=1) + 
    ggtitle('% of district-level cases in trial population')  + scale_color_manual(values=group.colors) + 
    scale_y_continuous(labels = formatC, limits=c(0,100)) 
## subs <- pf[, swctPT=='remPD_SF' & immunoDelay==21 & (trial == 'SWCT' & mod=='relabCoxME')]
## p.tmp <- p.tmp + geom_line(aes(vaccEff, caseTot), data = pf[subs], linetype = 4, col = 'orange', size = 1.3)
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Figure S9 - Cases by trial',typ), p.tmp, w = 9, h = 3.5)

subs <- pf[,vaccEff==0 & (delayUnit==0 | (trial=='SWCT')) & mod =='CoxME']
pf[subs, list(caseTot,caseC,caseV,pit,trial,delayUnit,vaccEff)]

####################################################################################################
## Get abstract & discussion #'s.
####################################################################################################
pf <- arrange(pf, immunoDelay, trial, propInTrial, mod)

subs <- pf[,mainAn==T & immunoDelay==21 & (trial=='SWCT' & mod=='relabCoxME') | (trial=='RCT' & ord=='TU' & mod=='CoxME')]
pf[subs & vaccEff==.9, list(trial, ord, delayUnit, mod, vaccGoodNAR, propInTrial,vaccEff, immunoDelay)]

subs <- pf[,immunoDelay==21 & mainAn==T  & mod %in% c('CoxME') & vaccEff>.7 &  (trial=='SWCT' | ( trial=='FRCT' & ord=='TU'))]
pf[subs, list(trial, ord, mod, vaccGoodNAR,cvr, biasNAR, propInTrial,vaccEff)]

####################################################################################################
## Table S1 - Bias & coverage
####################################################################################################
subs <- pf[,immunoDelay==21 & mainAn==T  & vaccEff>.7 & 
           ((mod %in% c('CoxME','GLMFclus','bootCoxME','bootGLMFclus','relabCoxME') & trial=='SWCT' ) | 
           (mod %in% c('CoxME') & trial=='RCT' & ord=='TU' ))]
pf[subs, list(trial, ord, mod, vaccGoodNAR,cvr, biasNAR, propInTrial,vaccEff)]
## note GLMFclus has better coverage here & less bias, but since Cox is more flexible for non-monotonic trends, let's stick with it
subs <- pf[,immunoDelay==21 & mainAn==T  & vaccEff>.7 & 
           ((mod %in% c('CoxME','relabCoxME') & trial=='SWCT' ) | 
           (mod %in% c('CoxME') & trial=='RCT' & ord=='TU' ))]
pf[subs, list(trial, ord, mod, vaccGoodNAR,cvr, biasNAR, propInTrial,vaccEff)]

tab1 <- pf[subs, list(trial, mod, pit, cvr, biasNAR, vaccEff)]
tab1[, list(trial, cvr, biasNAR)]
tab1[trial=='SWCT', cvr:= cvr[!is.na(cvr)], list(pit,vaccEff)]
tab1 <- tab1[!(trial=='SWCT' & mod=='CoxME')]
tab1[, c('cvr','biasNAR') := list(signif(cvr,4), signif(biasNAR,4))]
tab1 <- data.table(tab1[trial=='RCT',list(pit, cvr,biasNAR)], tab1[trial=='SWCT',list(cvr,biasNAR)])
tab1
write.csv(tab1, file='Results/biasCoverage.csv')

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
