if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

thing <- 'initDateSens'
batchdirnm <- file.path('BigResults',thing)
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

dparms <- c('trial','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit','immunoDelay','trialStartDate')
nbatch <- length(fls)
finInfoList <- finModList <- stopList <- parmsList <- list(NULL)
for(ii in 1:nbatch) {
    if(ii%%100 ==0) print(ii)
    ff <- fls[ii]
    if(exists('sim')) rm(sim)
    load(ff)
    if(exists('sim')) {
        sim$parms[['trialStartDate']] <- as.character(sim$parms[['trialStartDate']])
        parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
         tmpMod <- data.frame(nbatch = ii, sim$sim$finMods)
        finModList[[ii]] <- merge(tmpMod, sim$sim$finInfo, by = 'sim')
    }
}

parmsDT <- rbindlist(parmsList)
finTrials <- merge(rbindlist(finModList), parmsDT, by = c('nbatch'))
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]

finTrials[, sum(is.na(p)), mod]
finTrials[, sum(is.na(lci)), mod]
finTrials[, length(lci), list(propInTrial, mod)]
finTrials[mod=='coxME' & is.na(p), err:=1] ## sometimes cox returns NaNs, or partial NA's for certain values
finTrials$vaccEff <- as.numeric(finTrials$vaccEff)

## Determine if stopped
finTrials[grepl('boot',mod), stopped := lci > 0 | uci < 0]
finTrials[grepl('relab',mod), stopped := p < .025]
finTrials[!grepl('boot',mod) & !grepl('relab',mod), stopped := p < .05]
finTrials[, vaccGood := stopped==T &  mean > 0]
finTrials[, vaccBad := stopped==T &  mean < 0]
## Coverage
finTrials[, cvr := lci < vaccEff & uci > vaccEff]
finTrials[is.na(cvr), cvr := F]
finTrials[grepl('relab',mod), cvr := NA] # no CI's for perm test
## Bias
finTrials[, bias := mean - vaccEff]
## Reorder columns
front <- c('mod','stopped','vaccGood','vaccBad')
setcolorder(finTrials, c(front, setdiff(names(finTrials), front)))
back <- c('nbatch','sim')
setcolorder(finTrials, c(setdiff(names(finTrials), back), back))
save(finTrials, file=file.path('BigResults', paste0(thing, '.Rdata')))

load(file=file.path('BigResults',paste0(thing, '.Rdata')))

powFin <- summarise(group_by(finTrials, vaccEff, trial, propInTrial, ord, delayUnit, mod, immunoDelay,trialStartDate)
                    , nsim = length(stopped)
                    , stopped = mean(stopped)
                    , vaccGood = mean(vaccGood)
                    , cvr = mean(cvr)
                    , cvrNAR = mean(cvr, na.rm=T)
                    , bias = mean(bias)
                    , biasNAR = mean(bias, na.rm=T)
                    , vaccBad = mean(vaccBad)
                    , stoppedNAR = mean(stopped,na.rm=T)
                    , vaccGoodNAR = mean(vaccGood,na.rm=T)
                    , vaccBadNAR = mean(vaccBad,na.rm=T)
                    , meanErr = mean(err)
                    , meanBump = mean(bump)
                    , stopDay = mean(stopDay)
                    , caseTot = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    , caseC = mean(caseCXimmGrpEnd)
                    , caseV = mean(caseVXimmGrpEnd)
                    )
powFin[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]
powFin[,delayUnit:= as.numeric(levels(delayUnit)[delayUnit])]
powFin[,trial:=factor(trial)]
front <- c('mod','vaccEff','stopped','stoppedNAR','vaccGood','vaccGoodNAR','cvr','cvrNAR','bias','biasNAR',
'nsim','meanErr','propInTrial','vaccBad')
setcolorder(powFin, c(front, setdiff(names(powFin), front)))
pf <- data.table(powFin)
pf <- pf[!(trial=='FRCT' & delayUnit==0) & !(ord=='TU' & delayUnit==0)] ## redundant

pf$trialStartDate <- as.Date(pf$trialStartDate)
save(pf, file=file.path('Results',paste0('powFin_',thing,'.Rdata')))

## to delete a range of jobs
## qdel echo `seq -f "%.0f" 2282389 2282404`
