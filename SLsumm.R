if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSims')
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

dparms <- c('trial','ord','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize')
nbatch <- length(fls)
finList <- stopList <- parmsList <- list(NULL)
length(stopList) <- length(finList) <- length(parmsList) <- nbatch
for(ii in 1:nbatch) {
    if(ii%%100 ==0) print(ii)
    ff <- fls[ii]
    load(ff)
    parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
    stopList[[ii]] <- data.frame(nbatch = ii, sim$sim$stopPoints)
    finList[[ii]] <- data.frame(nbatch = ii, sim$sim$endFinRes)
}
parmsDT <- rbindlist(parmsList)
stopTrials <- merge(rbindlist(stopList), parmsDT, by = 'nbatch')
finTrials <- merge(rbindlist(finList), parmsDT, by = 'nbatch')
stopTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]
save(stopTrials, finTrials, file=file.path('Results','SLSumm.Rdata'))

load(file=file.path('Results','SLSumm.Rdata'))
 
powFin <- summarise(group_by(finTrials[sdLogIndiv==1], vaccEff, trial, propInTrial)
                    , stopped = mean(stopped)
                    , vaccGood = mean(stopped & vaccGood)
                    , stopDay = mean(stopDay)
                    , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC_stopActive = mean(caseCXimmGrpEnd)
                    , caseV_stopActive = mean(caseVXimmGrpEnd)
                    )

powStop <- summarise(group_by(stopTrials[sdLogIndiv==1], vaccEff, trial, propInTrial)
                    , stopped = mean(stopped)
                    , vaccGood = mean(stopped & vaccGood)
                    , stopDay = mean(stopDay)
                    , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC_stopActive = mean(caseCXimmGrpEnd)
                    , caseV_stopActive = mean(caseVXimmGrpEnd)
                    , totCase_finActive = mean(caseCXrandFinA + caseVXrandFinA)
                    , caseC_finActive = mean(caseCXrandFinA)          
                    , caseV_finActive = mean(caseVXrandFinA)
                    , totCase_fin = mean(caseCXrandFin + caseVXrandFin)          
                    , caseC_fin = mean(caseCXrandFin)
                    , caseV_fin = mean(caseVXrandFin)
                    )

powStop[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]
powFin[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]

## Power by efficacy & weekly decay rate, panels by weeklydecayvar
cols <- rainbow(4)
pdf('Figures/power by efficacy in SL propInTrial.pdf', w = 8, h = 6)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of subnational\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1), bty = 'n', main = main, las = 1)
    abline(h=.025)
    powFin[propInTrial==pit & vaccEff <=.85,
           lines(vaccEff, vaccGood, col = cols[as.numeric(trial)]), 
           by = list(trial)]
}
plot.new()
legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
title(main='24 week power', outer = T)
mtext('probability of rejecting the null hypothesis', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

## Power by efficacy & weekly decay rate, panels by weeklydecayvar
cols <- rainbow(4)
pdf('Figures/number cases by efficacy in SL propInTrial.pdf', w = 8, h = 6)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]
for(ii in 1:length(pits)) {
    pit <- pits[ii]
    main <- paste0('proportion of subnational\n cases in trial = ', signif(pit,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,130), bty = 'n', main = main, las = 1)
    abline(h=.025)
    powFin[propInTrial==pit & vaccEff <=.85,
           lines(vaccEff, totCase_stopActive, col = cols[as.numeric(trial)]), 
           by = list(trial)]
}
plot.new()
legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
title(main='24 week power', outer = T)
mtext('# cases', 1, 0, outer = T)
mtext('vaccine efficacy', 2, 0, outer = T)
graphics.off()

## Power by efficacy & weekly decay rate, panels by cluster level variation
cols <- rainbow(4)
pdf('Figures/power by efficacy wdr cvclus.pdf', w = 8, h = 6)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
vcls <- powFin[,unique(cvClus)]
for(ii in 1:length(vcls)) {
    vcl <- vcls[ii]
    main <- paste0('cluster-level coef of var = ', signif(vcl,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1), bty = 'n', main = main, las = 1)
    abline(h=.025)
    powFin[weeklyDecayVar==makeParms()$weeklyDecayVar & cvClus==vcl, lines(vaccEff, vaccGood, col = cols[as.numeric(trial)], lty = as.numeric(weeklyDecay)),
           by = list(cvClus, weeklyDecay, weeklyDecayVar,trial)]
    if(ii==3)    legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
    if(ii==4)    legend('topleft', leg=powFin[,levels(weeklyDecay)], lty = 1:powFin[,nlevels(weeklyDecay)], lwd = 2, bty = 'n', title = 'weekly decay rate')
}
title(main='35 week power', outer = T)
mtext('power', 1, 0, outer = T)
mtext('vaccine efficacy', 2, 0, outer = T)
graphics.off()


p1 <- simTrial(makeParms('RCT',small=F, clusSize=1, numClus = 20, weeklyDecay= .95, weeklyDecayVar = wdvs[2]))
nm <- 'default variation in EVD hazard \nfor trial participants'
plotHazT(p1, flnm='hazard trajectories' , main=nm )

