if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer)
##load(file=file.path('BigResults','powFin.Rdata'))

thing <- 'SLSims5'
load(file=file.path('BigResults',paste0('powFin_',thing,'.Rdata')))

head(pF[,1:8,with=F],50)

maxVE <- 1
maxPwr <- 1
cols <- c('black','orange','green')#,'blue')
## cols <- rainbow(4)
cols <- data.table(tri = pF[,levels(trial)], col = cols)
modtypes <- pF[,levels(mod)]

pdf('Figures/power SL propInTrial TU vs ord.pdf', w = 8, h = 6) ##, units = 'in', res = 200)
qplot(vaccEff, stoppedNAR, colour data = pF)
dev.off()


## power by order
##jpeg('Figures/power SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
pdf('Figures/power SL propInTrial TU vs ord.pdf', w = 8, h = 6) ##, units = 'in', res = 200)
for(mm in 1:length(modtypes)) {
    modTmp <- modtypes[mm]
    subs <- pF[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                   ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) &
                   vaccEff <=maxVE]
    powTmp <- pF[subs]
    powTmp[,trial:=factor(trial)]
    par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
    for(ii in 1:length(pits)) {
        pit <- pits[ii]
        main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
        plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,maxPwr), bty = 'n', main = main, las = 1)
        powTmp[propInTrial==pit, {
            lty <- which(c('none','TU')==ord)
            if(delayUnit==0) lty <- 3
            lines(vaccEff, vaccGoodNAR, col = cols[tri==trial, col], lty = lty)
        },
               by = list(trial, ord, delayUnit > 0)]
    }
    plot.new()
    tnmsTmp <- powTmp[,unique(trial)]
    legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
    legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
    title(main=paste0('24 week power; ',modTmp) , outer = T)
    mtext('probability of rejecting the null hypothesis', 2, 0, outer = T)
    mtext('vaccine efficacy', 1, 0, outer = T)
}
graphics.off()

## Type 1 Error
for(yvar in c('stopped','vaccGood','vaccBad')) {
    pdf(paste0('Figures/FalsePos', yvar, ' SL propInTrial TU vs ord.pdf'), w = 8, h = 6) ##, units = 'in', res = 200)
    ##jpeg('Figures/FalsePos SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
    par(lwd=2, mar = c(5,5,3,.5), mfrow = c(2,3), oma = c(1,1,1,0))
    for(mm in 1:length(modtypes)) {
        modTmp <- modtypes[mm]
        subs <- pF[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                       ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
                       & vaccEff ==0 ]
        powTmp <- pF[subs]
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
    pF[mod==modtypes[1] & propInTrial==pit, {
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
subs <- pF[, mod==modtypes[1] & !(delayUnit==0 & ord!='none') & 
               ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
               & vaccEff ==0 ]
powTmp <- pF[subs]
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
    pF[mod==modtypes[1] & propInTrial==pit, {
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
        subs <- pF[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                       ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
                       & vaccEff ==0 ]
        powTmp <- pF[subs]
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
