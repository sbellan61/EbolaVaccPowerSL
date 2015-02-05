if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)


dparms <- c('trial','ord','mu','varClus','sdLogIndiv','vaccEff','weeklyDecay','weeklyDecayVar')
nbatch <- length(fls)
finList <- stopList <- parmsList <- list(NULL)
length(stopList) <- length(finList) <- length(parmsList) <- nbatch
for(ii in 1:nbatch) {
    ff <- fls[ii]
    load(ff)
    parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
    stopList[[ii]] <- data.frame(nbatch = ii, sim$sim$stopPoints)
    finList[[ii]] <- data.frame(nbatch = ii, sim$sim$endFinRes)
}
parmsDT <- rbindlist(parmsList)
stopTrials <- merge(rbindlist(stopList), parmsDT, by = 'nbatch')
finTrials <- merge(rbindlist(finList), parmsDT, by = 'nbatch')
save(stopTrials, finTrials, file=file.path('Results','fpSumm.Rdata'))

finTrials[vaccEff==0 & sdLogIndiv==1 ,list(fp = mean(stopped)), list(varClus, weeklyDecay, sdLogIndiv, vaccEff, trial)]
finTrials[vaccEff==0 & sdLogIndiv==1 ,list(fp = mean(stopped & vaccGood)), list(varClus, weeklyDecay, sdLogIndiv, vaccEff, trial)]
stopTrials[vaccEff==0 & sdLogIndiv==1 ,list(fp = mean(stopped & vaccGood)), list(varClus, weeklyDecay, sdLogIndiv, vaccEff, trial)]

stopSummary <- summarise(group_by(stopTrials[stopped & vaccGood & sdLogIndiv==1], varClus, weeklyDecay, vaccEff, trial)
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
stopSummary <- data.table(stopSummary)
head(stopSummary)
stopSummary[trial=='RCT']


select(stopSummary, sdLogIndiv, vaccEff, trial, stopDay, totCase_stopActive, totCase_finActive, totCase_fin)
stopSummary[, caseC_fin]`

