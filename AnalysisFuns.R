## Construct survival data from waiting times
makeSurvDat <- function(parms,  whichDo='pop', browse=F) within(parms, {
    if(browse) browser()
    popTmp <- get(whichDo)
    popTmp$immuneDayThink <- popTmp[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    ## pre-immunity table
    stPre <- copy(popTmp) # st = survival table
    stPre$startDay <- 0
    stPre$endDay <- stPre[, pmin(immuneDayThink, infectDay)]
    stPre$infected <- stPre[ ,as.numeric(infectDay <= immuneDayThink)]
    stPre$immuneGrp <-  0     ## immuneGrp is variable used for analysis, not omnietient knowledge of vaccination/immune status
    stPre <- stPre[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, 
                         immuneDayThink, startDay, endDay, infectDay, infected, immuneGrp)]
    ## post-immunity table
    stPost <- copy(popTmp)[infectDay > immuneDayThink,]
    stPost$startDay <- stPost[,immuneDayThink]
    stPost$endDay   <-  stPost[,infectDay]
    stPost$infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
    stPost$immuneGrp <- 1
    stPost <- stPost[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, 
                           immuneDayThink, startDay, endDay, infectDay,  infected, immuneGrp)]
    nmSt <- paste0(sub('pop','',whichDo),'st') ## makes EVpop into EVst for example
    assign(nmSt, rbind(stPre, stPost)) ## combine tables
    rm(stPre, stPost, popTmp, nmSt)
}) ## careful with modifying parms, st depends on analysis a bit too (immunoDelayThink), so we can have different st for same popH

## Select subset of survival table to analyze
activeFXN <- function(parms, whichDo='st', browse=F) within(parms, { 
    if(browse) browser()
    ## for SWCT or unmatched CRCT always include all clusters in analysis because unvaccinated
    ## clusters are still considered to provide useful information from baseline
    stA <- copy(get(whichDo))
    stA$firstActive <- 0
    if(!includeAllControlPT) { ## remove person-time observed prior to post-refractory period from data
        if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
            stA[, firstActive := min(immuneDayThink), by = pair]
        if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
            stA[, firstActive := min(immuneDayThink), by = cluster]
    }
    stA <- stA[!endDay <= firstActive] ## remove inactive observation intervals
    stA[startDay < firstActive, startDay := firstActive] ## set accumulation of person time as when the cluster/pair is active
    nmStA <- paste0(sub('st','',whichDo),'stActive') ## makes EVpop into EVst for example
    assign(nmStA, stA)
    rm(stA, nmStA)
})
## p1 <- simTrial(makeParms('RCT', ord='BL', small=F), br=F)
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## s1$st[idByClus%in%1:2, list(indiv, cluster, pair, idByClus,immuneDayThink, startDay,endDay)]

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = parms$maxInfectDay+parms$hazIntUnit, whichDo = 'stActive') with(parms, {
    stTmp <- copy(get(whichDo))
    intervalNotStarted <- stTmp[,startDay] > censorDay
    stTmp <- stTmp[!intervalNotStarted,] 
    noInfectionBeforeCensor <- stTmp[,endDay] > censorDay
    stTmp[noInfectionBeforeCensor, infected:=0]
    stTmp[noInfectionBeforeCensor, endDay:=censorDay]
    stTmp[,perstime := (endDay-startDay)]
    stTmp <- stTmp[perstime > 0,] 
    return(stTmp)
})

summTrial <- function(st) list(summarise(group_by(st, cluster), sum(infected))
                               , summarise(group_by(st, cluster, immuneGrp), sum(infected))
                               , summarise(group_by(st, immuneGrp), sum(infected))
                               )

compileStopInfo <- function(minDay, vaccEffEst, tmp) {
    out <- data.table(stopDay=minDay
                      , mean = vaccEffEst['mean'], lci = vaccEffEst['lci'], uci = vaccEffEst['uci'], p = vaccEffEst['p']
                      , caseVaccEndTXImmuneGrp = tmp[immuneGrp==1, sum(infected)]
                      , caseContEndTXImmuneGrp = tmp[immuneGrp==0, sum(infected)]
                      ## , ptVaccEndTXImmuneGrp = tmp[immuneGrp==1, sum(perstime)]
                      ## , ptContEndTXImmuneGrp = tmp[immuneGrp==0, sum(perstime)]
                      )
    ## out$hazVaccEndTXImmuneGrp <- out[, caseVaccEndTXImmuneGrp / ptVaccEndTXImmuneGrp / yearToDays]
    ## out$hazContEndTXImmuneGrp <- out[, caseContEndTXImmuneGrp / ptContEndTXImmuneGrp / yearToDays]
    ## out$ptRatioEndTXImmuneGrp <- out[, ptContEndTXImmuneGrp / ptVaccEndTXImmuneGrp]
    out$stopped <- out[, p<.05 & !is.na(p)]
    out$vaccGood <- NA
    out[stopped==T, vaccGood:= lci > 0]
    out <- setcolorder(out, c("stopped", "vaccGood", "stopDay", "mean", "lci", "uci", "p"
                              , "caseVaccEndTXImmuneGrp", "caseContEndTXImmuneGrp"
                              ## , "ptVaccEndTXImmuneGrp", "ptContEndTXImmuneGrp"
                              ## , "hazVaccEndTXImmuneGrp", "hazContEndTXImmuneGrp",  "ptRatioEndTXImmuneGrp"
                              ))
    out <- as.data.frame(out)
    return(out)
}

## Check whether stopping point has been reached at intervals
seqStop <- function(parms, start = parms$immunoDelayThink + 14, checkIncrement = 7, verbose = 0, maxDay = parms$maxInfectDay) {
    trialOngoing <- T
    checkDay <- start
    first <- T
    while(trialOngoing) {
        if(verbose>1) browser()
        if(verbose>0) print(checkDay)
        tmp <- censSurvDat(parms, checkDay)
        vaccEffEst <- try(doCoxPH(tmp), silent=T) ## converting midDay to days from months
        ## if cox model has enough info to converge check for stopping criteria
        if(!inherits(vaccEffEst, 'try-error') & !is.nan(vaccEffEst['p'])) { 
            newout <- compileStopInfo(checkDay, vaccEffEst, tmp) 
            if(first) out <- newout else out <- rbind(out, newout)
            first <- F
            if(!is.na(newout['p']))
                if(newout['p'] < .05)
                    trialOngoing <- F
            if(checkDay > maxDay) trialOngoing <- F
        }
        checkDay <- checkDay + 7
    }
    rownames(out) <- NULL
    out <- as.data.table(out)
    parms$weeklyAns <- out
    parms$endTrialDay <- tail(out$stopDay,1)
    parms$vaccEffEst <- vaccEffEst
    return(parms)
}

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, check=F, returnAll = F, verbose=0) {
    set.seed(seed)
    casesXVaccRandGrpList <- casesXPT_ImmuneList <- weeklyAnsList <- list()
    length(casesXVaccRandGrpList) <- length(casesXPT_ImmuneList) <- length(weeklyAnsList) <- N
    stopPoints <- data.frame(NULL)
    for(ii in 1:N) {
        if(verbose>1) browser()
        res <- simTrial(parms)
        res <- makeSurvDat(res)
        res <- activeFXN(res)
        res <- seqStop(res)
        browser()
        res <- endT(res,T)
        res <- makeCaseSummary(res)

        ## active cases at end of trial
        ## total cases at final
        ## active cases at final
        stopPt <- as.data.frame(tail(res$weeklyAns,1))
        stopPt <- c(stopPt
                    , contCasesFinalActive = res$casesXVaccRandGrp[type=='EVstActive', contCases] 
                    , vaccCasesFinalActive = res$casesXVaccRandGrp[type=='EVstActive', vaccCases]
                    , contCasesFinal = res$casesXVaccRandGrp[type=='EVst', contCases] 
                    , vaccCasesFinal = res$casesXVaccRandGrp[type=='EVst', vaccCases]
                    )
        stopPoints <- rbind(stopPoints, stopPt)
        if(returnAll) {
            weeklyAnsList[[ii]] <- as.data.frame(res$weeklyAns)
            casesXVaccRandGrpList[[ii]] <- as.data.frame(res$casesXVaccRandGrp)
            casesXPT_ImmuneList[[ii]] <- as.data.frame(res$casesXPT_Immune)
        }
    }
    rownames(stopPoints) <- NULL
    if(returnAll)
        return(list(stopPoints, weeklyAnsList, casesXVaccRandGrpList, casesXPT_ImmuneList))
    if(!returnAll)
        return(stopPoints)
}
