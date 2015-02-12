if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSims3')
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
#fls <- list.files(batchdirnm, pattern='SL-2', full.names = T)
length(fls)

dparms <- c('trial','ord','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize','ord','reordLag')
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

## load(file=file.path('Results','SLSumm.Rdata'))
 
powFin <- summarise(group_by(finTrials[sdLogIndiv==1], vaccEff, trial, propInTrial, ord)
                    , nsim = length(stopped)
                    , stopped = mean(stopped)
                    , vaccGood = mean(stopped & vaccGood)
                    , stopDay = mean(stopDay)
                    , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC_stopActive = mean(caseCXimmGrpEnd)
                    , caseV_stopActive = mean(caseVXimmGrpEnd)
                    )


powFin[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]

    powFin[vaccEff <=.85 & (trial=='SWCT' | (trial %in% c('FRCT','RCT') & ord=='TU'))]

powFin <- powFin[trial!='CRCT']
powFin[, trial:= factor(trial, levels = unique(trial))]

## Power by efficacy & weekly decay rate, panels by weeklydecayvar
cols <- rainbow(4)
##pdf('Figures/power by efficacy in SL propInTrial.pdf', w = 8, h = 6)
jpeg('Figures/power by efficacy in SL propInTrial.jpg', w = 8, h = 6, units = 'in', res = 200)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,.6), bty = 'n', main = main, las = 1)
    abline(h=.025, lty = 3)
    powFin[propInTrial==pit & vaccEff <=.85 & (trial=='SWCT' | (trial %in% c('FRCT','RCT') & ord=='TU')),
           lines(vaccEff, vaccGood, col = cols[as.numeric(trial)]), 
           by = list(trial)]
}
plot.new()
legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
title(main='24 week power', outer = T)
mtext('probability of rejecting the null hypothesis', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

cols <- rainbow(4)
pdf('Figures/number cases by efficacy in SL propInTrial.pdf', w = 8, h = 6)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of district-level\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,130), bty = 'n', main = main, las = 1)
    abline(h=.025)
    powFin[propInTrial==pit & vaccEff <=.85,
           lines(vaccEff, totCase_stopActive, col = cols[as.numeric(trial)]), 
           by = list(trial)]
}
plot.new()
legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
title(main='24 week power', outer = T)
mtext('# cases in trial', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()



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
