if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(boot)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

thing <- 'FalsePosFluct'
batchdirnm <- file.path('BigResults',thing)
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'Phenom', weeklyDecay = .9, cvWeeklyDecay = .5, cvClus = 1.5, cvClusTime = 0.5, numClus = 20))

dparms <- c('trial','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit','immunoDelay','trialStartDate'
            , 'weeklyDecay', 'cvWeeklyDecay', 'cvClus', 'cvClusTime'
            )
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

## Simulations with less than 10 cases are considered to not have any power
finTrials$tooSmall <- finTrials[, (caseCXimmGrpEnd + caseVXimmGrpEnd) < 10]
finTrials[tooSmall==T, c('vaccGood','vaccBad','stopped') := F]
finTrials[tooSmall==T, c('lci','uci','p') := list(-Inf,1,1)]
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
## Bias, must be done on RH/(RH+1) scale to deal with Inf & 0 RH's
finTrials$RH <- finTrials[, 1-mean]
finTrials$PHU <- finTrials[, RH/(RH+1)] ## proportion of hazard unavoidable even with vaccination
finTrials[RH==Inf, PHU:=1] ## otheriwse gives NaN for Inf/(Inf+1)
finTrials[,list(vaccEff,mean,PHU)] ## NEED TO AVERAGE BIAS ON PHU scale 

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
                    , mean = mean(mean)
                    , meanNAR = mean(mean, na.rm=T)
                    , vaccBad = mean(vaccBad)
                    , stoppedNAR = mean(stopped,na.rm=T)
                    , vaccGoodNAR = mean(vaccGood,na.rm=T)
                    , vaccBadNAR = mean(vaccBad,na.rm=T)
                    , PHUNAR = mean(PHU, na.rm=T)
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
## Get bias from means done on a PHU scale
powFin$RHxPHUNAR <- powFin[, -PHUNAR/(PHUNAR-1)]
powFin$meanXPHUNAR <- powFin[, 1 - RHxPHUNAR]
powFin$biasNAR <- powFin[, meanXPHUNAR - vaccEff]
powFin[vaccEff==.7, list(meanNAR,meanXPHUNAR,vaccEff, biasNAR)]

## Formatting stuff
front <- c('mod','vaccEff','stoppedNAR','vaccGoodNAR','cvrNAR','biasNAR',
'nsim','meanErr','propInTrial','vaccBad','cvr','stopped','vaccGood')
setcolorder(powFin, c(front, setdiff(names(powFin), front)))
pf <- data.table(powFin)
pf <- pf[!(trial=='FRCT' & delayUnit==0) & !(ord=='TU' & delayUnit==0)] ## redundant
pf$trialStartDate <- as.Date(pf$trialStartDate)
pf[mod=='coxME', mod:='CoxME']
pf$mod <- factor(pf$mod, levels=unique(pf$mod))
pf$order <- pf$ord
pf[delayUnit==0, order:='simultaneous instant']
pf$design <- pf$trial
levels(pf$design)[levels(pf$design) == 'SWCT'] <- 'SWT'
levels(pf$order)[2] <- 'time-updated'
pf[, immunoDelay:=as.numeric(levels(immunoDelay)[immunoDelay])]
pf[, pit:=factor(paste0(propInTrial*100,'%'))]
pf[, pit:=factor(pit, levels = c('2.5%','5%','7.5%','10%'), ordered = T)]
baseMods <- c('Cox PH Frailty'
              , 'Poisson GLM\n no cluster effects'
              , 'Poisson GLM \nwith fixed effects by cluster')
pf$model <- pf$mod
levels(pf$model) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))

save(pf, file=file.path('Results',paste0('powFin_',thing,'.Rdata')))

## to delete a range of jobs
## qdel echo `seq -f "%.0f" 2282389 2282404`


