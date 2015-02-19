if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)
##load(file=file.path('BigResults','powFin.Rdata'))

thing <- 'SLSims5'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))
pf <- pF
pf <- pf[!(trial=='FRCT' & delayUnit==0) & !(ord=='TU' & delayUnit==0)]
pf[,vaccEff:=as.numeric(vaccEff)]
pf[delayUnit==0, ord:='simultaneous instant']
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.03 & mod=='coxME']
setnames(pf,c('ord','trial'),c('order','design'))
levels(pf$order)[2] <- 'time-updated'
pf[, pit:=factor(paste0(propInTrial*100,'%'), levels = c('3%','5%','10%'), ordered=T)]

baseMods <- c('Cox PH Frailty'
              , 'Poisson GLM\n no cluster effects'
              , 'Poisson GLM \nwith fixed effects by cluster')

pf$modLb <- pf$mod
levels(pf$modLb) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))

pdf('Figures/power SL propInTrial TU vs ord.pdf', w = 8, h = 5) ##, units = 'in', res = 200)
for(modTmp in levels(pf$modLb)) {
    thax <- element_text(colour = 'black', size = 8)
    p.tmp <- ggplot(pf[modLb==modTmp], aes(vaccEff, vaccGoodNAR, colour=design, linetype=order)) + 
        scale_x_continuous(labels = formatC, limits=c(0,1)) +  
        ##scale_y_continuous(labels = formatC, limits=c(0,1), minor_breaks = seq(0,1,.05)) +
        scale_y_log10(labels = formatC, limits=c(0.01,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8)) +
            theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
                  axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5)) +
                      xlab('vaccine efficacy') + ylab('probability of rejecting null') + 
                          scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                              geom_line(size=1) + facet_wrap(~propInTrial) + 
                                  ggtitle('power by expected % of district-level cases in trial population') +
                        geom_hline(yintercept=.025, color='black', linetype='dotted')
    print(p.tmp)
    grid.text(modTmp,x=unit(.85,"npc"),y=unit(0.85,"npc"), gp=gpar(fontsize=10))
}
graphics.off()


pdf('Figures/ SL propInTrial TU vs ord.pdf', w = 8, h = 5) ##, units = 'in', res = 200)
for(modTmp in levels(pf$modLb)) {
    thax <- element_text(colour = 'black', size = 8)
    p.tmp <- ggplot(pf[mod==modTmp], aes(vaccEff, vaccGoodNAR, colour=design, linetype=order)) + 
        scale_x_continuous(labels = formatC, limits=c(0,1)) +  
        ##scale_y_continuous(labels = formatC, limits=c(0,1), minor_breaks = seq(0,1,.05)) +
        scale_y_log10(labels = formatC, limits=c(0.01,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8)) +
            theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
                  axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
                  panel.margin = unit(2, "lines")) +
                      xlab('vaccine efficacy') + ylab('probability of rejecting null') + 
                          scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                              geom_line(size=1) + facet_wrap(~propInTrial) + 
                                  ggtitle('power by expected % of district-level cases in trial population') +
                        geom_hline(yintercept=.025, color='black', linetype='dotted')
    print(p.tmp)
    grid.text(modTmp,x=unit(.85,"npc"),y=unit(0.85,"npc"), gp=gpar(fontsize=10))
}
graphics.off()

## Type 1 Error
for(yvar in c('stopped','vaccGood','vaccBad')) {
    pdf(paste0('Figures/FalsePos', yvar, ' SL propInTrial TU vs ord.pdf'), w = 8, h = 6) ##, units = 'in', res = 200)
    ##jpeg('Figures/FalsePos SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
    par(lwd=2, mar = c(5,5,3,.5), mfrow = c(2,3), oma = c(1,1,1,0))
    for(mm in 1:length(modtypes)) {
        modTmp <- modtypes[mm]
        subs <- pf[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                       ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
                       & vaccEff ==0 ]
        powTmp <- pf[subs]
        powTmp[,trial:=factor(trial)]
        plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0.03,.1), ylim = c(0,.2), bty = 'n', main = '', las = 1, xaxt='n')
        axis(1, at = c(0,.03, .05, .1))
        powTmp[,
               {
                   lty <- which(c('none','TU')==ord)
                   if(delayUnit==0) lty <- 3
                   lines(propInTrial, get(yvar), col = cols[tri==trial, col], lty = lty, type = 'b')
               },
               by = list(trial, ord, delayUnit > 0)]
        abline(h=ifelse(yvar=='stopped',.05,.025), lty = 3)
        title(main=paste0(modTmp), outer = F, line = 0)
    }
    plot.new()
    legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
    legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
    mtext('Type I error rate', 2, -1, outer = T)
    mtext(yvar, 3, -1, outer = T)
    mtext('proportion of district-level cases in trial population', 1, -1, outer = T)
    graphics.off()
}

## cases by order
maxCases <- 120
jpeg('Figures/cases SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,maxCases), bty = 'n', main = main, las = 1)
    pf[mod==modtypes[1] & propInTrial==pit, {
        lty <- which(c('none','TU')==ord)
        if(delayUnit==0) lty <- 3
        lines(vaccEff, totCase_stopActive, col = cols[tri==trial, col], lty = lty)
    },
           by = list(trial, ord, delayUnit > 0)]
}
plot.new()
legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
title(main='24 week power', outer = T)
mtext('# cases in trial', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

####################################################################################################
## numsims
pdf(paste0('Figures/NSIMS SL fp.pdf'), w = 8, h = 6) ##, units = 'in', res = 200)
par(lwd=2, mar = c(5,5,3,.5), oma = c(1,1,1,0))
layout(matrix(1:2,1,2), w = c(3,1))
subs <- pf[, mod==modtypes[1] & !(delayUnit==0 & ord!='none') & 
               ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
               & vaccEff ==0 ]
powTmp <- pf[subs]
powTmp[,trial:=factor(trial)]
plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0.03,.1), ylim = c(0,1200), bty = 'n', main = '', las = 1, xaxt='n')
axis(1, at = c(0,.03, .05, .1))
powTmp[,
       {
           lty <- which(c('none','TU')==ord)
           if(delayUnit==0) lty <- 3
           lines(propInTrial, nsim, col = cols[tri==trial, col], lty = lty, type = 'b')
       },
       by = list(trial, ord, delayUnit > 0)]
par(mar=rep(0,4))
plot.new()
legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
mtext('# simulations', 2, -1, outer = T)
mtext('proportion of district-level cases in trial population', 1, -1, outer = T)
graphics.off()

jpeg('Figures/nsims SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1200), bty = 'n', main = main, las = 1)
    pf[mod==modtypes[1] & propInTrial==pit, {
        lty <- which(c('none','TU')==ord)
        if(delayUnit==0) lty <- 3
        lines(vaccEff, nsim, col = cols[tri==trial, col], lty = lty)
    },
           by = list(trial, ord, delayUnit > 0)]
}
plot.new()
legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
title(main='24 week power', outer = T)
mtext('# cases in trial', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

####################################################################################################
## SWCT Type 1 Error
for(yvar in c('stopped','vaccGood','vaccBad')) {
    pdf(paste0('Figures/FalsePos', yvar, thing, ' SL propInTrial TU vs ord.pdf'), w = 8, h = 7) ##, units = 'in', res = 200)
    ##jpeg('Figures/FalsePos SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
    par(lwd=2, mar = c(5,5,3,.5), mfrow = c(3,3), oma = c(1,1,1,0))
    for(mm in 1:length(modtypes)) {
        modTmp <- modtypes[mm]
        subs <- pf[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                       ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
                       & vaccEff ==0 ]
        powTmp <- pf[subs]
        powTmp[,trial:=factor(trial)]
        plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0.03,.1), ylim = c(0,.2), bty = 'n', main = '', las = 1, xaxt='n')
        axis(1, at = c(.03, .05, .1))
        powTmp[,
               {
                   lty <- which(c('none','TU')==ord)
                   if(delayUnit==0) lty <- 3
                   lines(propInTrial, get(yvar), lty = lty, type = 'b')
               },
               by = list(trial, ord, delayUnit > 0)]
        abline(h=ifelse(yvar=='stopped',.05,.025), lty = 3)
        title(main=paste0(modTmp), outer = F, line = 0)
    }
    mtext('Type I error rate', 2, -1, outer = T)
    mtext(yvar, 3, 0, outer = T)
    mtext('proportion of district-level cases in trial population', 1, -1, outer = T)
    graphics.off()
}
