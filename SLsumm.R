if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSimsFP')
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
finTrials <- merge(rbindlist(finList), parmsDT, by = c('nbatch'))
## stopTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[mod=='bootCoxME', stopped := lci > 0 | uci < 0]
finTrials[mod=='bootCoxME', vaccGood := lci > 0]
finTrials[mod=='bootCoxME' & !stopped, vaccGood := NA]
finTrials[mod=='relabCoxME', stopped := p < .05]
finTrials[mod=='relabCoxME', vaccGood := p < .05 & mean > 0 ]
save(finTrials, file=file.path('BigResults','SLSummFP.Rdata'))

load(file=file.path('Results','SLSummFP.Rdata'))

powFin <- summarise(group_by(finTrials[sdLogIndiv==1], vaccEff, trial, propInTrial, ord, delayUnit, mod)
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
powFin[,vaccEff:= as.numeric(levels(vaccEff)[vaccEff])]
powFin[,trial:=factor(trial)]
powFin <- powFin[!(mod=='GEEClusAR1' & trial %in% c('RCT','FRCT'))]
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]



maxVE <- 1
maxPwr <- 1
cols <- c('black','brown','green','blue')
cols <- rainbow(4)
cols <- data.table(tri = powFin[,levels(trial)], col = cols)
modtypes <- powFin[,levels(mod)]

## power by order
##jpeg('Figures/power SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
pdf('Figures/power SL propInTrial TU vs ord.pdf', w = 8, h = 6) ##, units = 'in', res = 200)
for(mm in 1:length(modtypes)) {
    modTmp <- modtypes[mm]
    subs <- powFin[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                   ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) &
                   vaccEff <=maxVE]
    powTmp <- powFin[subs]
    powTmp[,trial:=factor(trial)]
    par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
    for(ii in 1:length(pits)) {
        pit <- pits[ii]
        main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
        plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,maxPwr), bty = 'n', main = main, las = 1)
        powTmp[propInTrial==pit, {
            lty <- which(c('none','TU')==ord)
            if(delayUnit==0) lty <- 3
            lines(vaccEff, vaccGood, col = cols[tri==trial, col], lty = lty)
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
pdf('Figures/FalsePos SL propInTrial TU vs ord.pdf', w = 8, h = 6) ##, units = 'in', res = 200)
##jpeg('Figures/FalsePos SL propInTrial TU vs ord.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mar = c(5,5,3,.5), mfrow = c(2,3), oma = c(1,1,0,0))
for(mm in 1:length(modtypes)) {
    modTmp <- modtypes[mm]
    subs <- powFin[, mod==modTmp & !(delayUnit==0 & ord!='none') & 
                   ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
                   & vaccEff ==0 ]
    powTmp <- powFin[subs]
    powTmp[,trial:=factor(trial)]
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0.03,.1), ylim = c(0,.2), bty = 'n', main = '', las = 1, xaxt='n')
    axis(1, at = c(0,.03, .05, .1))
    ##    abline(h=.025, lty = 3)
    powTmp[,
           {
               lty <- which(c('none','TU')==ord)
               if(delayUnit==0) lty <- 3
               lines(propInTrial, vaccGood, col = cols[tri==trial, col], lty = lty, type = 'b')
           },
           by = list(trial, ord, delayUnit > 0)]
    abline(h=.025, lty = 3)
    title(main=paste0(modTmp), outer = F, line = 0)
}
plot.new()
legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
mtext('Type I error rate', 2, -1, outer = T)
mtext('proportion of district-level cases in trial population', 1, -1, outer = T)
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
    powTmp[mod==modtypes[1] & propInTrial==pit, {
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

## to delete a range of jobs
## qdel echo `seq -f "%.0f" 2282389 2282404`
