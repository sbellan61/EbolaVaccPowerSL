if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSimsFPbump')
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

dparms <- c('trial','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit')
nbatch <- length(fls)
finList <- stopList <- parmsList <- list(NULL)
length(stopList) <- length(finList) <- length(parmsList) <- nbatch
for(ii in 1:nbatch) {
    if(ii%%100 ==0) print(ii)
    ff <- fls[ii]
    if(exists('sim')) rm(sim)
    load(ff)
    if(exists('sim')) {
        parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
        finList[[ii]] <- data.frame(nbatch = ii, sim$sim$finPoint)
    }
}

parmsDT <- rbindlist(parmsList)
## stopTrials <- merge(rbindlist(stopList), parmsDT, by = 'nbatch')
finTrials <- merge(rbindlist(finList), parmsDT, by = c('nbatch'))
## stopTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]
finTrials[mod=='bootCoxME', stopped := lci > 0 | uci < 0]
finTrials[mod=='bootCoxME' & stopped==T, vaccGood := lci > 0]
finTrials[mod=='bootCoxME' & !stopped, vaccGood := NA]
finTrials[mod=='relabCoxME', stopped := p < .05]
finTrials[mod=='relabCoxME' & stopped==T, vaccGood := p < .05 & mean > 0 ]
finTrials[mod=='relabCoxME' & !stopped, vaccGood := NA]
finTrials[is.na(stopped), stopped := F]
finTrials[is.na(vaccGood), stopped := F]
finTrials$vaccBad <- finTrials[, !vaccGood]
finTrials[is.na(vaccGood), vaccGood := F] ## fill in Falses for unstopped trials
finTrials[is.na(vaccBad), vaccBad := F]
save(finTrials, file=file.path('BigResults','SLSummFP.Rdata'))

load(file=file.path('BigResults','SLSummFP.Rdata'))

powFin <- summarise(group_by(finTrials[sdLogIndiv==1], vaccEff, trial, propInTrial, ord, delayUnit, mod)
                    , nsim = length(stopped)
                    , stopped = mean(stopped)
                    , vaccGood = mean(vaccGood)
                    , vaccBad = mean(vaccBad)
                    , stopDay = mean(stopDay)
                    , totCase_stopActive = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC_stopActive = mean(caseCXimmGrpEnd)
                    , caseV_stopActive = mean(caseVXimmGrpEnd)
                    )
powFin[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]
powFin[,delayUnit:= as.numeric(levels(delayUnit)[delayUnit])]
##powFin[,vaccEff:= as.numeric(levels(vaccEff)[vaccEff])]
powFin[,vaccEff:= as.numeric(vaccEff)]
powFin[,trial:=factor(trial)]
powFin <- powFin[!(mod=='GEEClusAR1' & trial %in% c('RCT','FRCT'))]
pits <- powFin[,unique(propInTrial)]
pits <- pits[order(pits)]
save(powFin, file=file.path('BigResults','powFin.Rdata'))

## to delete a range of jobs
## qdel echo `seq -f "%.0f" 2282389 2282404`

batchdirnm <- file.path('BigResults','SLSimsFPbump','Routs')
routs <- list.files(batchdirnm, pattern='.Rout', full.names = T)
find.error <- function(fl, msg = 'Error') sum(grepl(msg, readLines(fl)))>0
errs <- unlist(lapply(routs, find.error))
rerr <- routs[errs]
head(rerr)

errs <- unlist(lapply(rerr, find.error), msg = 'Error in if (effMean > 0) pval <- 1 - pval')
