if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)
##load(file=file.path('BigResults','powFin.Rdata'))
percent <- function(x) paste0(formatC(x*100), '%')
labs <- c('','log')

thing <- 'SLSimsFinal'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.025 & mod=='CoxME']
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
group.colors <- c(RCT = "#333BFF", FRCT = "#CC6600", SWT ="#9633FF")
group.colors[c(1,3,2)] <- gg_color_hue(3)
group.colors['SWCT'] <- 'orange'
pf$trial <- factor(pf$trial, levels=levels(pf$trial)[c(2,1,3)])
pf[, biasNAR:=biasNAR/vaccEff]
levels(pf$order)[1] <- 'random'
levels(pf$order)[2] <- 'risk-prioritized'


####################################################################################################
## them for ms
thax <- element_text(colour = 'black', size = 8)
thsb <- theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
              axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
              axis.line = element_line(), axis.ticks = element_line(color='black'),
              panel.margin = unit(1, "lines"), legend.key.height=unit(2,"line")
              , strip.background = element_rect(fill = NA)
              ,legend.position = 'right'
              , axis.line = element_blank()
              ,panel.grid.major = element_blank()
              , panel.grid.minor = element_blank()
              ,panel.border = element_blank()
              ,panel.background = element_blank()
              , legend.background =  element_blank()
              , legend.key =  element_blank()
              , legend.key.width=unit(2,"line")
              ## ,legend.justification=c(1,0), legend.position=c(1,0)
              )
theme_set(theme_grey(base_size = 12))
##thsb <- thsb + theme_bw()#

pf$modLab2 <- pf$mod
levels(pf$modLab2)[levels(pf$modLab2) %in% c('CoxME','bootCoxME','relabCoxME')] <- c('(A) CoxPH', '(B) Bootstrap', '(C) Permutation')

####################################################################################################
## Figure 2 - Type I errors
for(jj in 1:2) {
    subs <- pf[,trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
    subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
    p.tmp <- ggplot(pf[subs], 
                    aes(propInTrial, stoppedNAR, colour=trial, linetype=order)) + thsb +
                        scale_x_continuous(labels = percent, limits=c(.025,.1), minor_breaks=NULL, breaks = c(.025,.05,.075,.1)) +  
                            xlab('% of district-level cases in trial population') + ylab('False Positive Rate') + 
                                scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                    geom_hline(yintercept=.05, color='dark gray', size = 1) +
                                        geom_line(size=1) + facet_wrap(~modLab2, scales = "free_y") + scale_color_manual(values=group.colors)
    if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.1))
    if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1), limits = c(.005,.1), minor_breaks=NULL)
    ggsave(paste0('Figures/Fig 2A -',labs[jj],'Type I SL.png'), p.tmp, w = 8.5, h = 3.5)
    ggsave(paste0('Figures/Fig 2A -',labs[jj],'Type I SL.pdf'), p.tmp, w = 8.5, h = 3.5)
}

####################################################################################################
## Figure 2B - Type I errors
for(jj in 1:2) {
thax <- element_text(colour = 'black', size = 8)
subs <- pf[,trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=model, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.025,.1), minor_breaks=NULL, breaks = c(.025,.05,.075,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~trial, scales = "free_y") + 
        scale_color_manual(values=group.colors)
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25), limits = c(.005,.25), minor_breaks=NULL) 
ggsave(paste0('Figures/Fig 2B -',labs[jj],'Type I SL.png'), p.tmp, w = 6.5, h = 3.5)
ggsave(paste0('Figures/Fig 2B -',labs[jj],'Type I SL.pdf'), p.tmp, w = 6.5, h = 3.5)
}

####################################################################################################
## Figure 4 - Power
subs <- pf[, immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT','FRCT') & mod =='CoxME'))]
for(jj in 1:2) {
    thax <- element_text(colour = 'black', size = 8)
    p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=trial, linetype=order)) + thsb +
        scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
            xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
                scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                    geom_line(size=1) + facet_wrap(~pit, scales = "free_y",nrow=1) + 
                        ggtitle('expected % of district-level cases in trial population') + 
                            scale_color_manual(values=group.colors)
    if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
    if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, limits=c(0.01,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8))
    ggsave(paste0('Figures/Fig 4 -',labs[jj],'Power SL.png'), p.tmp, w = 9, h = 4)
    ggsave(paste0('Figures/Fig 4 -',labs[jj],'Power SL.pdf'), p.tmp, w = 9, h = 4)
}

####################################################################################################
## Figure 5 - Power by immunedelay
subs <- pf[, propInTrial==.05 & vaccEff == .9 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial=='RCT' & order=='risk-prioritized' & mod =='CoxME'))]
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[subs], aes(immunoDelay, vaccGoodNAR, colour=trial, linetype=order)) + thsb +
         #   scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
                          xlab('delay until protective efficacy') + ylab('power to detect effective vaccine') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
scale_color_manual(values=group.colors) +
                                  geom_line(size=1)  #facet_wrap(~immunoDelay, scales = "free_y",nrow=1) + 
 p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
ggsave(paste0('Figures/Fig 5 - Power by seroconversion delay SL.png'), p.tmp, w = 5, h = 4)
ggsave(paste0('Figures/Fig 5 - Power by seroconversion delay SL.pdf'), p.tmp, w = 5, h = 4)



####################################################################################################
## Figure SX - Power by model
pf$modelLab <- pf$model
pf$class <- 'parametric'
pf[grepl('boot',mod), class:='bootstrap']
pf[grepl('relab',mod), class:='permutation']
pf$class <- factor(pf$class, levels = c("parametric", "bootstrap", "permutation"))
pf$model <- gsub('boot','',pf$mod)
pf$model <- factor(gsub('relab','',pf$model))
levels(pf$model)[1:2] <- c('CoxPH','Poisson regression')
subs <- pf[, propInTrial==.05  & immunoDelay==21 & (trial == 'SWCT' | (trial=='RCT' & order=='risk-prioritized' & class!='bootstrap'))]
subs <- subs & pf[,model %in% c('CoxPH','Poisson regression')]
thax <- element_text(colour = 'black', size = 8)
for(jj in 1:2) {
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=class, linetype=model)) +
    thsb +
    scale_x_continuous(labels = formatC, limits=c(0,1), breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('probability of rejecting null \nthat vaccine is not effective') + 
    scale_linetype_manual(values=1:3) + ## scale_color_manual(values=group.colors) +
    geom_line(size=.9)  + facet_wrap(~trial, scales = "free_y") #,nrow=1) + 
#p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25, .5, 1), limits = c(.005,1), minor_breaks=NULL) 
ggsave(paste0('Figures/Fig SX', labs[jj], ' - Power by model delay SL.png'), p.tmp, w = 6, h = 4)
ggsave(paste0('Figures/Fig SX', labs[jj], ' - Power by model delay SL.pdf'), p.tmp, w = 6, h = 4)
}

##################################################
## Cases
subs <- pf[, immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT','FRCT') & mod =='CoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, caseTot, colour=trial, linetype=order)) + thsb + 
    scale_x_continuous(labels = formatC, limits=c(0,.9), breaks = pf[,unique(vaccEff)],minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('# EVD cases') +
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit, scale='free_y', nrow=1) + 
    ggtitle('% of district-level cases in trial population')  + scale_color_manual(values=group.colors)
p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,120)) #
print(p.tmp)
ggsave(paste0('Figures/Cases SL ALL.png'), p.tmp, w = 9, h = 3.5)
ggsave(paste0('Figures/Cases SL ALL.pdf'), p.tmp, w = 9, h = 3.5)

subs <- pf[,vaccEff==0 & (delayUnit==0 | (trial=='SWCT')) & mod =='CoxME']
pf[subs, list(caseTot,caseC,caseV,pit,trial,delayUnit,vaccEff)]

####################################################################################################
## abstract #'s and tables
pf <- arrange(pf, immunoDelay, trial, propInTrial, mod)

pf[immunoDelay==21 & mod %in% c('relabCoxME') & vaccEff>.8 &  ((trial=='SWCT' & ord=='random')| ( trial=='RCT' & ord=='TU')),
list(trial, ord, delayUnit, mod, vaccGoodNAR, propInTrial,vaccEff)]

pf[immunoDelay==21 & mod %in% c('CoxME') & vaccEff==0 &  ((trial=='SWCT' & ord=='random')| ( trial=='RCT' & ord=='TU')),
list(trial, ord, delayUnit, mod, stoppedNAR, propInTrial,vaccEff)]

pf[mod %in% c('CoxME') & vaccEff>.5 &  ((trial=='SWCT' & ord=='random')| ( trial=='FRCT' & ord=='TU')),
list(trial, ord, mod, vaccGoodNAR,cvr, propInTrial,vaccEff)]

## Chosen models


subs <- pf[, vaccEff %in% c(0,.9) & immunoDelay==21 & #propInTrial == .05 & 
           ((mod %in% c('CoxME','relabCoxME') &  (trial=='SWCT' & ord=='random') ) |
           (mod %in% c('CoxME') & trial=='RCT' & ord=='TU' ))]
tab1 <- pf[subs, list(trial, mod, pit, cvr, biasNAR, vaccEff)]
tab1
tab1[, list(trial, cvr, biasNAR)]
tab1[trial=='SWCT', cvr:= cvr[!is.na(cvr)], list(pit,vaccEff)]
tab1 <- tab1[!(trial=='SWCT' & mod=='CoxME')]
tab1[, c('cvr','biasNAR') := list(signif(cvr,2), signif(biasNAR,2))]
tab1

##dcast.data.table(tab1, trial + mod + pit ~ vaccEff, value.var=c('cvr','biasNAR'))
tab1 <- data.table(tab1[trial=='RCT',list(pit, cvr,biasNAR)], tab1[trial=='SWCT',list(cvr,biasNAR)])
tab1
write.csv(tab1, file='Results/biasCoverage.csv')




#################################################
## Power
for(jj in 1:2) {
    pdf(paste0('Figures/',labs[jj], 'power SL ALL.pdf'), w = 8, h = 4) ##, units = 'in', res = 200)
    for(modTmp in levels(pf$model)) {
        subs <- pf[, model==modTmp & immunoDelay==21 & !(trial%in%c('FRCT','RCT') & grepl('boot',mod))]
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=trial, linetype=order)) + 
            scale_x_continuous(labels = formatC, limits=c(0,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
                theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
                      axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
                      panel.margin = unit(2, "lines"))  + scale_color_manual(values=group.colors) + 
                          xlab('vaccine efficacy') + ylab('probability of rejecting null') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                  geom_line(size=1) + facet_wrap(~pit, nrow=1) + 
                                      ggtitle('power by expected % of district-level cases in trial population') +
                                          geom_hline(yintercept=.025, color='black', linetype='dotted')
        if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), minor_breaks = seq(0,1,.05))
        if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, limits=c(0.0005,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8))
        print(p.tmp)
        grid.text(modTmp,x=unit(.85,"npc"),y=unit(0.85,"npc"), gp=gpar(fontsize=10))
    }
    graphics.off()
}

##################################################
## Type I
labs <- c('','log')
jj <- 1

for(jj in 1:2) {
pdf(paste0('Figures/',labs[jj],'Type I SL.pdf'), w = 8, h = 8) ##, units = 'in', res = 200)
thax <- element_text(colour = 'black', size = 8)
subs <- pf[,trial %in% c('SWCT','RCT') & vaccEff==0 & immunoDelay==21 & !(trial=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], aes(propInTrial, stoppedNAR, colour=trial, linetype=order)) + 
    scale_x_continuous(labels = percent, limits=c(.025,.1), breaks = c(.025,.05,.075,.1)) +  
##    scale_y_continuous(labels = formatC, limits=c(0,.3)) + #, minor_breaks = seq(0,1,.05)) +
    theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
          panel.margin = unit(2, "lines")) + scale_color_manual(values=group.colors) +
    xlab('% of distrinct-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~model) + 
geom_hline(yintercept=.05, color='black', linetype='dotted', size = 1.2)
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC)#, limits=c(0,.3))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .2, .3, .4))#,  limits=c(.01,.3))
print(p.tmp)
graphics.off()
}
