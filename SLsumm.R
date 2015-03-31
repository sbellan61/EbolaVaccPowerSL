####################################################################################################
## Collect results of simulation ensembles from HPC cluster and put into analyzeable data tables.
####################################################################################################
## Code base accompanying:
## 
## Bellan, SE, JRC Pulliam, CAB Pearson, DChampredon, SJ Fox, L Skrip, AP Galvani, M Gambhir, BA
## Lopman, TC Porco, LA Meyers, J Dushoff (2015). The statistical power and validity of Ebola
## vaccine trials in Sierra Leone: A simulation study of trial design and analysis. _Lancet
## Infectious Diseases_.
##
## Steve Bellan, March 2015
## License at bottom.
####################################################################################################

sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

thing <- 'All'
batchdirnm <- file.path('BigResults',thing)
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
fls <- fls[grepl('pit', fls)]
length(fls)

dparms <- c('trial','sdLogIndiv','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit','immunoDelay','trialStartDate'
            , 'weeklyDecay', 'cvWeeklyDecay', 'cvClus', 'cvClusTime', 'remStartFin', 'remProtDel'
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

parmsDT <- rbindlist(parmsList, use.names = T, fill = T)

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

powFin <- summarise(group_by(finTrials, vaccEff, trial, propInTrial, ord, delayUnit, mod, immunoDelay,trialStartDate
                             , weeklyDecay, cvWeeklyDecay, cvClus, remProtDel, remStartFin) #, cvClusTime)
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

####################################################################################################
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/sbellan61/EbolaVaccPowerSL
####################################################################################################
