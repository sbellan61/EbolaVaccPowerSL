if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','FalsePosSims2')
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

dparms <- c('trial','ord','mu','hazSL', 'cvClus','sdLogIndiv','vaccEff','weeklyDecay','cvWeeklyDecay','delayUnit')
nbatch <- length(fls)
finList <- stopList <- parmsList <- list(NULL)
length(stopList) <- length(finList) <- length(parmsList) <- nbatch
for(ii in 1:nbatch) {
    if(ii%%100 ==0) print(ii)
    ff <- fls[ii]
    load(ff)
    parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
    ## stopList[[ii]] <- data.frame(nbatch = ii, sim$sim$stopPoints)
    finList[[ii]] <- data.frame(nbatch = ii, sim$sim$finPoint)
}
parmsDT <- rbindlist(parmsList)
## stopTrials <- merge(rbindlist(stopList), parmsDT, by = 'nbatch')
finTrials <- merge(rbindlist(finList), parmsDT, by = 'nbatch')
## stopTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]
## stopTrials[,cvClus := sqrt(as.numeric(levels(varClus)[varClus])) / as.numeric(levels(mu)[mu])]
save(finTrials, file=file.path('Results','fpSumm.Rdata'))

load(file=file.path('Results','fpSumm.Rdata'))
 
powFin <- summarise(group_by(finTrials[sdLogIndiv==1], cvClus, weeklyDecay, cvWeeklyDecay, vaccEff, trial,  ord)
                    , nsim = length(stopped)
                    , stopped = mean(stopped)
                    , vaccGood = mean(stopped & vaccGood)
                    , stopDay = mean(stopDay)
                    , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC_stopActive = mean(caseCXimmGrpEnd)
                    , caseV_stopActive = mean(caseVXimmGrpEnd)
                    )
powFin[,cvWeeklyDecay:= as.numeric(levels(cvWeeklyDecay)[cvWeeklyDecay])]
powFin[,weeklyDecay:= as.numeric(levels(weeklyDecay)[weeklyDecay])]
powFin[,cvClus:= as.numeric(levels(cvClus)[cvClus])]
wds <- powFin[,unique(weeklyDecay)]
wds <- wds[order(wds)]
cvcls <- powFin[,unique(cvClus)]
cvcls <- cvcls[order(cvcls)]
cvwds <- powFin[,unique(cvWeeklyDecay)]
cvwds <- cvwds[order(cvwds)]

maxVE <- 1
maxPwr <- 1

## Power by efficacy & weekly decay rate, panels by weeklydecayvar
cols <- rainbow(4)
jpeg('Figures/power by efficacy Exphaz wdr.jpg', w = 8, h = 6, units = 'in', res = 300)
par(lwd=2, mfrow=c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
for(ii in 1:length(cvcls)) {
    cvwd <- cvwds[ii]
    main <- paste0('weekly decay CV = ', signif(cvwd,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1), bty = 'n', main = main, las = 1)
    abline(h=.025)
    powFin[cvWeeklyDecay==cvwd & vaccEff <= .85,
           lines(vaccEff, vaccGood, col = cols[as.numeric(trial)], lty = which(cvClus == cvcls)),
           by = list(cvClus, weeklyDecay, cvWeeklyDecay,trial)]
}
plot.new()
legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
legend('left', leg=cvcls, lty = 1:length(cvcls), lwd = 2, bty = 'n', 
       title = 'cv of baseline cluster hazard')
title(main='24 week power', outer = T)
mtext('probability of rejecting the null', 2, 0, outer = T)
mtext('vaccine efficacy', 1, 0, outer = T)
graphics.off()

powFin[vaccEff==0]

## ## Power by efficacy & weekly decay rate, panels by weeklydecayvar STOP TIME
## cols <- rainbow(4)
## pdf('Figures/power at STOP TIME by efficacy wdr.pdf', w = 8, h = 6)
## par(lwd=2, mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
## cvwds <- powStop[,unique(cvWeeklyDecay)]
## for(ii in 1) {
##     cvwd <- cvwds[ii]
##     cvwdst <- as.numeric(levels(cvwd)[cvwd])
##     main <- paste0('weekly decay variance = ', signif(cvwdst,3))
##     plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1), bty = 'n', main = main, las = 1)
##     abline(h=.025)
##     powStop[cvWeeklyDecay==cvwd & signif(cvClus,2)==signif(with(makeParms(), sqrt(cvClus)/mu),2),
##            lines(vaccEff, vaccGood, col = cols[as.numeric(trial)], lty = as.numeric(weeklyDecay)),
##            by = list(cvClus, weeklyDecay, cvWeeklyDecay,trial)]
##         legend('topleft', leg=powStop[,levels(trial)], col = cols, lwd = 2, bty = 'n')
##         legend('left', leg=powStop[,levels(weeklyDecay)], lty = 1:powStop[,nlevels(weeklyDecay)], lwd = 2, bty = 'n', title = 'weekly decay rate')
## }
## title(main='power checking weekly', outer = T)
## mtext('power', 1, 0, outer = T)
## mtext('vaccine efficacy', 2, 0, outer = T)
## graphics.off()


## Power by efficacy & weekly decay rate, panels by weeklydecayvar
cols <- rainbow(4)
pdf('Figures/power by efficacy wdr wdrvar.pdf', w = 8, h = 6)
par(lwd=2, mfrow = c(2,2), mar = c(3,3,3,.5), oma = c(1.5,1.5,1.5,0))
cvwds <- powFin[,unique(cvWeeklyDecay)]
cvwds <- cvwds[order(cvwds)]
wds <- unique(powFin[,weeklyDecay])
wds <- wds[rev(order(wds))][-4]
for(ii in 1:length(cvwds)) {
    cvwd <- cvwds[ii]
    main <- paste0('weekly decay variance = ', signif(cvwd,3))
    plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1), bty = 'n', main = main, las = 1)
    abline(h=.025)
    powFin[trial!='FRCT' & weeklyDecay > .9 & cvWeeklyDecay==cvwd & signif(cvClus,2)==signif(with(makeParms(), sqrt(cvClus)/mu),2),
           lines(vaccEff, vaccGood, col = cols[as.numeric(trial)], lty = which(weeklyDecay==wds)),
           by = list(cvClus, weeklyDecay, cvWeeklyDecay,trial)]
    if(ii==3)    legend('topleft', leg=powFin[,levels(trial)][-3], col = cols[-3], lwd = 2, bty = 'n')
    if(ii==1)    legend('topleft', leg=wds, lty = 1:length(wds), lwd = 2, bty = 'n', title = 'weekly decay rate')
}
title(main='35 week power', outer = T)
mtext('power', 1, 0, outer = T)
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
    powFin[cvWeeklyDecay==makeParms()$cvWeeklyDecay & cvClus==vcl, lines(vaccEff, vaccGood, col = cols[as.numeric(trial)], lty = as.numeric(weeklyDecay)),
           by = list(cvClus, weeklyDecay, cvWeeklyDecay,trial)]
    if(ii==3)    legend('topleft', leg=powFin[,levels(trial)], col = cols, lwd = 2, bty = 'n')
    if(ii==4)    legend('topleft', leg=powFin[,levels(weeklyDecay)], lty = 1:powFin[,nlevels(weeklyDecay)], lwd = 2, bty = 'n', title = 'weekly decay rate')
}
title(main='35 week power', outer = T)
mtext('power', 1, 0, outer = T)
mtext('vaccine efficacy', 2, 0, outer = T)
graphics.off()


p1 <- simTrial(makeParms('RCT',small=F, clusSize=1, numClus = 20, weeklyDecay= .95, cvWeeklyDecay = cvwds[2]))
nm <- 'default variation in EVD hazard \nfor trial participants'
plotHazT(p1, flnm='hazard trajectories' , main=nm )


## finTrials[vaccEff==0 & sdLogIndiv==1 ,list(fp = mean(stopped)), list(cvClus, weeklyDecay, sdLogIndiv, vaccEff, trial)]
## finTrials[vaccEff==0 & sdLogIndiv==1 ,list(fp = mean(stopped & vaccGood)), list(cvClus, weeklyDecay, sdLogIndiv, vaccEff, trial)]
## stopTrials[vaccEff==0 & sdLogIndiv==1 ,list(fp = mean(stopped & vaccGood)), list(cvClus, weeklyDecay, sdLogIndiv, vaccEff, trial)]

## stopSumam
## ry <- summarise(group_by(stopTrials[stopped & vaccGood & sdLogIndiv==1], cvClus, weeklyDecay, cvWeeklyDecay, vaccEff, trial)
##                          , stopDay = mean(stopDay)
##                          , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
##                          , caseC_stopActive = mean(caseCXimmGrpEnd)
##                          , caseV_stopActive = mean(caseVXimmGrpEnd)
##                          , totCase_finActive = mean(caseCXrandFinA + caseVXrandFinA)
##                          , caseC_finActive = mean(caseCXrandFinA)          
##                          , caseV_finActive = mean(caseVXrandFinA)
##                          , totCase_fin = mean(caseCXrandFin + caseVXrandFin)          
##                          , caseC_fin = mean(caseCXrandFin)
##                          , caseV_fin = mean(caseVXrandFin)
##                          )

## stopSummary <- data.table(stopSummary)
## head(stopSummary)
## stopSummary[trial=='RCT']


## select(stopSummary, sdLogIndiv, vaccEff, trial, stopDay, totCase_stopActive, totCase_finActive, totCase_fin)
## stopSummary[, caseC_fin]`


## powStop <- summarise(group_by(stopTrials[sdLogIndiv==1], cvClus, weeklyDecay, cvWeeklyDecay, vaccEff, trial)
##                     , stopped = mean(stopped)
##                     , vaccGood = mean(stopped & vaccGood)
##                     , stopDay = mean(stopDay)
##                     , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
##                     , caseC_stopActive = mean(caseCXimmGrpEnd)
##                     , caseV_stopActive = mean(caseVXimmGrpEnd)
##                     ## , totCase_finActive = mean(caseCXrandFinA + caseVXrandFinA)
##                     ## , caseC_finActive = mean(caseCXrandFinA)          
##                     ## , caseV_finActive = mean(caseVXrandFinA)
##                     ## , totCase_fin = mean(caseCXrandFin + caseVXrandFin)          
##                     ## , caseC_fin = mean(caseCXrandFin)
##                     ## , caseV_fin = mean(caseVXrandFin)
##                     )
## powStop[, weeklyDecay:=factor(weeklyDecay, levels = c('1','0.98','0.95','0.9'))]

## powStop[,cvWeeklyDecay:= as.numeric(levels(cvWeeklyDecay)[cvWeeklyDecay])]
