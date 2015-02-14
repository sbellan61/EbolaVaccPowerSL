if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSims3')
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
#fls <- list.files(batchdirnm, pattern='SL-2', full.names = T)
length(fls)

dparms <- c('trial','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit')
nbatch <- length(fls)
finList <- stopList <- parmsList <- list(NULL)
length(stopList) <- length(finList) <- length(parmsList) <- nbatch
for(ii in 1:nbatch) {
    if(ii%%100 ==0) print(ii)
    ff <- fls[ii]
    load(ff)
    parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
##    stopList[[ii]] <- data.frame(nbatch = ii, sim$sim$stopPoints)
    finList[[ii]] <- data.frame(nbatch = ii, sim$sim$finPoint)
}
parmsDT <- rbindlist(parmsList)
## stopTrials <- merge(rbindlist(stopList), parmsDT, by = 'nbatch')
finTrials <- merge(rbindlist(finList), parmsDT, by = 'nbatch')
## stopTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]
save(finTrials, file=file.path('Results','SLSumm3.Rdata'))

load(file=file.path('Results','SLSumm3.Rdata'))

head(finTrials[lci==-Inf & trial=='RCT' & vaccEff==.9, nbatch])

finTrials[is.na(vaccGood) & trial=='RCT' & vaccEff==.9 & propInTrial==.1, list(trial, vaccGood, mean, lci, p, caseCXimmGrpEnd, caseVXimmGrpEnd)]
finTrials[is.na(vaccGood) & trial=='RCT' & vaccEff==.9 & propInTrial==.1, mean(caseCXimmGrpEnd)]

powFin <- summarise(group_by(finTrials[sdLogIndiv==1], vaccEff, trial, propInTrial, ord, delayUnit)
                    , nsim = length(stopped)
                    , stopped = mean(stopped)
                    , vaccGood = mean(stopped & vaccGood)
                    , stopDay = mean(stopDay)
                    , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC_stopActive = mean(caseCXimmGrpEnd)
                    , caseV_stopActive = mean(caseVXimmGrpEnd)
                    )
powFin[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]
powFin[,delayUnit:= as.numeric(levels(delayUnit)[delayUnit])]
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]

maxVE <- 1
maxPwr <- 1
cols <- c('black','brown','green','blue')
cols <- rainbow(3)

## power by order
subs <- powFin[, !(delayUnit==0 & ord!='none') & ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT')) & vaccEff <=maxVE ]
powTmp <- powFin[subs]
powTmp[,trial:=factor(trial)]
jpeg('Figures/power SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,maxPwr), bty = 'n', main = main, las = 1)
    ##    abline(h=.025, lty = 3)
    powTmp[propInTrial==pit, {
        lty <- which(c('none','TU')==ord)
        if(delayUnit==0) lty <- 3
        lines(vaccEff, vaccGood, col = cols[as.numeric(trial)], lty = lty)
    },
           by = list(trial, ord, delayUnit > 0)]
}
plot.new()
tnmsTmp <- powTmp[,unique(trial)]
legend('topleft', leg=tnmsTmp, col = cols, lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
title(main='24 week power', outer = T)
mtext('probability of rejecting the null hypothesis', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

## cases by order
maxCases <- 120
jpeg('Figures/cases SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,maxCases), bty = 'n', main = main, las = 1)
    ##    abline(h=.025, lty = 3)
    powTmp[propInTrial==pit, {
        lty <- which(c('none','TU')==ord)
        if(delayUnit==0) lty <- 3
        lines(vaccEff, totCase_stopActive, col = cols[as.numeric(trial)], lty = lty)
    },
           by = list(trial, ord, delayUnit > 0)]
}
plot.new()
legend('topleft', leg=powTmp[,levels(trial)], col = cols, lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
title(main='24 week power', outer = T)
mtext('# cases in trial', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

## Type 1 Error
subs <- powFin[, !(delayUnit==0 & ord!='none') & ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT')) & vaccEff ==0 ]
powTmp <- powFin[subs]
powTmp[,trial:=factor(trial)]
jpeg('Figures/FalsePos SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mar = c(5,5,1,.5))
plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0.03,.1), ylim = c(0,.07), bty = 'n', main = '', las = 1, xaxt='n')
axis(1, at = c(0,.03, .05, .1))
##    abline(h=.025, lty = 3)
powTmp[,
       {
           lty <- which(c('none','TU')==ord)
           if(delayUnit==0) lty <- 3
           lines(propInTrial, vaccGood, col = cols[as.numeric(trial)], lty = lty, type = 'b')
       },
       by = list(trial, ord, delayUnit > 0)]
abline(h=.025, lty = 3)
## plot.new()
legend('topleft', leg=powTmp[,levels(trial)], col = cols, lwd = 2, bty = 'n')
legend('topright', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
## title(main='24 week power', outer = T)
mtext('Type I error rate', 2, 4)#, outer = T)
mtext('proportion of district-level cases in trial population', 1, 3)#, outer = T)
graphics.off()
 

## 
## pdf('Figures/number cases by efficacy in SL propInTrial.pdf', w = 8, h = 6)
## par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
## pits <- powFin[,unique(propInTrial)]
## pits <- pits[order(pits)]
## for(ii in 1:length(pits)) {
##     pit <- pits[ii]
##     main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
##     plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,130), bty = 'n', main = main, las = 1)
##     abline(h=.025)
##     powFin[propInTrial==pit & vaccEff <=maxVE,
##            lines(vaccEff, totCase_stopActive, col = cols[as.numeric(trial)]), 
##            by = list(trial)]
## }
## plot.new()
## legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
## title(main='24 week power', outer = T)
## mtext('# cases in trial', 2, 0, outer = T)
## mtext('vaccine efficacy', 1, 0, outer = T)
## graphics.off()



## powStop <- summarise(group_by(stopTrials[sdLogIndiv==1], vaccEff, trial, propInTrial)
##                     , stopped = mean(stopped)
##                     , vaccGood = mean(stopped & vaccGood)
##                     , stopDay = mean(stopDay)
##                     , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
##                     , caseC_stopActive = mean(caseCXimmGrpEnd)
##                     , caseV_stopActive = mean(caseVXimmGrpEnd)
##                     , totCase_finActive = mean(caseCXrandFinA + caseVXrandFinA)
##                     , caseC_finActive = mean(caseCXrandFinA)          
##                     , caseV_finActive = mean(caseVXrandFinA)
##                     , totCase_fin = mean(caseCXrandFin + caseVXrandFin)          
##                     , caseC_fin = mean(caseCXrandFin)
##                     , caseV_fin = mean(caseVXrandFin)
##                     )
## powStop[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]



## ## Random, no vacc delays
## subs <- powFin[, trial %in% c('RCT','CRCT') & ord=='none' & delayUnit==0 & vaccEff <=maxVE ]
## 
## jpeg('Figures/power SL propInTrial no ord no delay.jpg', w = 8, h = 6, units = 'in', res = 200)
## par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
## for(ii in 1:length(pits)) {
##     pit <- pits[ii]
##     main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
##     plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,maxPwr), bty = 'n', main = main, las = 1)
## ##    abline(h=.025, lty = 3)
##     powFin[propInTrial==pit & subs,
##            lines(vaccEff, vaccGood, col = cols[as.numeric(trial)]), 
##            by = list(trial)]
## }
## plot.new()
## legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
## title(main='24 week power', outer = T)
## mtext('probability of rejecting the null hypothesis', 2, 0, outer = T)
## mtext('vaccine efficacy', 1, 0, outer = T)
## graphics.off()

## boostrap resampling procedure
## showing gee/glmm break
## noising up right panel 
