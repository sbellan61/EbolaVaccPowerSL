if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)
##load(file=file.path('BigResults','powFin.Rdata'))
percent <- function(x) paste0(formatC(x*100), '%')
labs <- c('','log')

thing <- 'SLSimsFinal'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))

pf[delayUnit==0, ord:='simultaneous instant']
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.03 & mod=='coxME']
setnames(pf,'ord','order')
pf$design <- pf$trial
levels(pf$design) <- c('RCT','FRCT','SWT')
levels(pf$order)[2] <- 'time-updated'
pf[, pit:=factor(paste0(propInTrial*100,'%'), levels = c('3%','5%','10%'), ordered=T)]

baseMods <- c('Cox PH Frailty'
              , 'Poisson GLM\n no cluster effects'
              , 'Poisson GLM \nwith fixed effects by cluster')

pf$model <- pf$mod
levels(pf$model) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))

##################################################
## Power
for(jj in 1:2) {
    pdf(paste0('Figures/',labs[jj], 'power SL ALL.pdf'), w = 8, h = 5) ##, units = 'in', res = 200)
    for(modTmp in levels(pf$model)) {
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[model==modTmp], aes(vaccEff, vaccGoodNAR, colour=design, linetype=order)) + 
            scale_x_continuous(labels = formatC, limits=c(0,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
                theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
                      axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
                      panel.margin = unit(2, "lines")) +
                          xlab('vaccine efficacy') + ylab('probability of rejecting null') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                  geom_line(size=1) + facet_wrap(~pit) + 
                                      ggtitle('power by expected % of district-level cases in trial population') +
                                          geom_hline(yintercept=.025, color='black', linetype='dotted')
        if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), minor_breaks = seq(0,1,.05))
        if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, limits=c(0.01,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8))
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
p.tmp <- ggplot(pf[design %in% c('SWCT','FRCT') & vaccEff==0], aes(propInTrial, stoppedNAR, colour=design, linetype=order)) + 
    scale_x_continuous(labels = percent, limits=c(.03,.1), breaks = c(.03,.05,.1)) +  
##    scale_y_continuous(labels = formatC, limits=c(0,.3)) + #, minor_breaks = seq(0,1,.05)) +
    theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
          panel.margin = unit(2, "lines")) +
    xlab('% of distrinct-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~model) + 
geom_hline(yintercept=.05, color='black', linetype='dotted', size = 1.2)
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC)#, limits=c(0,.3))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .2, .3, .4))#,  limits=c(.01,.3))
print(p.tmp)
graphics.off()
}



####################################################################################################
## them for ms
thsb <- theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
          axis.line = element_line(), axis.ticks = element_line(color='black'),
          panel.margin = unit(1, "lines"), legend.key.height=unit(2,"line")
,legend.position = 'right'
#              ,legend.justification=c(1,0), legend.position=c(1,0)
)
theme_set(theme_grey(base_size = 12))

####################################################################################################
## Figure 2 - Type I errors
for(jj in 1:2) {
thax <- element_text(colour = 'black', size = 8)
subs <- pf[,design %in% c('SWT','RCT') & vaccEff==0 & mod %in% c('coxME','bootCoxME','relabCoxME')]
subs <- subs & pf[,!(design=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=design, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.03,.1), minor_breaks=NULL, breaks = c(.03,.05,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~model, scales = "free_y")
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25), limits = c(.01,.25), minor_breaks=NULL) 
ggsave(paste0('Figures/Fig 2A -',labs[jj],'Type I SL.png'), p.tmp, w = 8.5, h = 3.5)
}


####################################################################################################
## Figure 2B - Type I errors
for(jj in 1:2) {
thax <- element_text(colour = 'black', size = 8)
subs <- pf[,design %in% c('SWT','RCT') & vaccEff==0 & mod %in% c('coxME','bootCoxME','relabCoxME')]
subs <- subs & pf[,!(design=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=model, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.03,.1), minor_breaks=NULL, breaks = c(.03,.05,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~design, scales = "free_y")
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25), limits = c(.01,.25), minor_breaks=NULL) 
ggsave(paste0('Figures/Fig 2B -',labs[jj],'Type I SL.png'), p.tmp, w = 6.5, h = 3.5)
}

####################################################################################################
## Figure 4 - Power
subs <- pf[, (design == 'SWT' & mod=='relabCoxME') | (design %in% c('RCT','FRCT') & mod =='coxME')]
for(jj in 1:2) {
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=design, linetype=order)) + thsb +
            scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
                          xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                  geom_line(size=1) + facet_wrap(~pit, scales = "free_y") + 
                                      ggtitle('expected % of district-level cases in trial population')
        if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
        if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, limits=c(0.01,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8))
ggsave(paste0('Figures/Fig 4 -',labs[jj],'Power SL.png'), p.tmp, w = 8, h = 4)
}


##################################################
## Cases
modTmp <- 'coxME'
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[mod==modTmp], aes(vaccEff, caseTot, colour=design, linetype=order)) + thsb + 
    scale_x_continuous(labels = formatC, limits=c(0,.9), breaks = pf[,unique(vaccEff)],minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('# EVD cases') +
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit, scale='free_y') + 
    ggtitle('% of district-level cases in trial population') +
    geom_hline(yintercept=.025, color='black', linetype='dotted')
p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,100)) #
print(p.tmp)
ggsave(paste0('Figures/Cases SL ALL.png'), p.tmp, w = 8, h = 3.5)

##################################################
## check # of simulations run, not too many crashed
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[mod==modTmp], aes(vaccEff, nsim, colour=design, linetype=order)) + 
    scale_x_continuous(labels = formatC, limits=c(0,1), breaks = pf[,unique(vaccEff)],minor_breaks=NULL) +  
    theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5)) +
    xlab('vaccine efficacy') + ylab('# simulations') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit) + 
    ggtitle('% of district-level cases in trial population') +
    geom_hline(yintercept=.025, color='black', linetype='dotted')
p.tmp <- p.tmp + scale_y_continuous(labels = formatC) #
print(p.tmp)
ggsave(paste0('Figures/numSims SL ALL.png'), p.tmp, w = 8, h = 3.5)

subs <- pf[,vaccEff==0 & (delayUnit==0 | (trial=='SWCT')) & mod =='coxME']
pf[subs, list(caseTot,caseC,caseV,pit,trial,delayUnit,vaccEff)]
