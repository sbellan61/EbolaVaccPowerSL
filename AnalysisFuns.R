## Construct survival data from waiting times
makeSurvDat <- function(parms,  whichDo='pop') within(parms, {
    if(verbose ==1.5) browser()
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
    nmSt <- paste0('st', sub('pop','',whichDo)) ## makes EVpop into EVst for example
    assign(nmSt, rbind(stPre, stPost)) ## combine tables
    rm(stPre, stPost, popTmp, nmSt)
}) ## careful with modifying parms, st depends on analysis a bit too (immunoDelayThink), so we can have different st for same popH

## Select subset of survival table to analyze
activeFXN <- function(parms, whichDo='st') within(parms, { 
    if(verbose ==1.6) browser()
    ## for SWCT or unmatched CRCT always include all clusters in analysis because unvaccinated
    ## clusters are still considered to provide useful information from baseline
    stA <- copy(get(whichDo))
    stA$firstActive <- 0
    if(!includeAllControlPT) { ## remove person-time observed prior to post-refractory period from data
        if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
            stA[, firstActive := min(immuneDayThink), by = pair]
        if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
            stA[, firstActive := min(immuneDayThink), by = cluster]
        if(trial=='SWCT') {## inactive during protective delay; active only when there exists both vacc & unvacc person-time observed
            firstDayAnyoneImmune <- stA[, min(immuneDayThink)]
            lastDayAnyoneNotImmune <- stA[, max(immuneDayThink)] - 1
            stA <- stA[!endDay <= firstDayAnyoneImmune] ## remove inactive observation intervals at beggining of trial
            stA <- stA[!startDay >= lastDayAnyoneNotImmune] ## remove inactive observation intervals at end of trial
            rm(firstDayAnyoneImmune, lastDayAnyoneNotImmune)
            ## remove person-time completely contained within protective delay (should only remove cluster 1's 0-immunedaythink person-time
            stA <- stA[!(startDay >= vaccDay & endDay <= immuneDayThink)]
            ## right-truncate at vaccine date person-time intervals that starts before vaccine date and ends
            ## after vaccine date (i.e. ignore person-time within protective delay)
            stA[startDay <= vaccDay & endDay >= vaccDay, endDay := vaccDay]
            ## stA[idByClus==1, list(cluster, vaccDay, immuneDayThink, startDay, endDay)] ## run to see how person-time is distributed between cluster
        }
    }
    stA <- stA[!endDay <= firstActive] ## remove inactive observation intervals (does nothing for SWCT)
    stA[startDay < firstActive, startDay := firstActive] ## set accumulation of person time as when the cluster/pair is active (does nothing for SWCT)
    nmStA <- paste0('stActive', sub('st','',whichDo)) ## makes EVpop into EVst for example
    assign(nmStA, stA)
    rm(stA, nmStA)
})
## p1 <- simTrial(makeParms('RCT', ord='BL', small=F), br=F)
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## s1$st[idByClus%in%1:2, list(indiv, cluster, pair, idByClus,immuneDayThink, startDay,endDay)]

## Restructure for GEE/GLMM with weekly observations of each cluster.
makeGEEDat <- function(parms, whichDo='popH') within(parms, {
    if(verbose ==1.7) browser()
    popHTmp <- get(whichDo)
    popHTmp$immuneDayThink <- popHTmp[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    popHTmp$infectDayRCens <- popHTmp$infectDay
    popHTmp[infectDay==Inf, infectDayRCens := NA]
    popHTmp$immuneGrp <- 0
    popHTmp[day >= immuneDayThink, immuneGrp := 1]
    popHTmp$firstActive <- 0
    if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
        popHTmp[, firstActive := min(immuneDayThink), by = pair]
    if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
        popHTmp[, firstActive := min(immuneDayThink), cluster]
    popHTmp$active <- popHTmp[,day>=firstActive]
    clusD <- popHTmp[active==T, list(cases = sum(!is.na(infectDayRCens))), list(cluster, day, immuneGrp)]
    clusD <- clusD[order(cluster, day)]
    clusD <- mutate(group_by(clusD, cluster), atRisk = clusSize - c(0, cumsum(cases[-length(cases)])))
    nmSt <- paste0('clusDat', sub('popH','',whichDo)) ## makes EVpop into EVst for example
    assign(nmSt, clusD) ## combine tables
    rm(popHTmp, clusD, nmSt)
})

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = parms$maxInfectDay+parms$hazIntUnit, whichDo = 'stActive') with(parms, {
    if(verbose==2.7) browser()
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

## Check whether stopping point has been reached at intervals
seqStop <- function(parms, start = parms$immunoDelayThink + 14, checkIncrement = 7, minCases = 15,
                    fullSeq = F, maxDay = parms$maxInfectDay) {
    trialOngoing <- T
    checkDay <- start
    first <- T
    while(trialOngoing) {
        if(parms$verbose>1) browser()
        if(parms$verbose>0) print(checkDay)
        tmp <- censSurvDat(parms, checkDay)
        vaccEffEst <- try(doCoxPH(tmp), silent=T) ## converting midDay to days from months
        ## if cox model has enough info to converge check for stopping criteria
        if(!inherits(vaccEffEst, 'try-error') & !is.nan(vaccEffEst['p'])) { 
            newout <- compileStopInfo(checkDay, vaccEffEst, tmp) 
            if(first) out <- newout else out <- rbind(out, newout)
            first <- F
            numCases <- newout['caseCXimmGrpEnd'] + newout['caseVXimmGrpEnd']
            if(!fullSeq & !is.na(newout['p']) & numCases > minCases)
                if(newout['p'] < .05)
                    trialOngoing <- F
        }
        checkDay <- checkDay + 7
        if(checkDay > maxDay) trialOngoing <- F
    }
    rownames(out) <- NULL
    out <- as.data.table(out)
    parms$weeklyAns <- out
    parms$endTrialDay <- tail(out$stopDay,1)
    parms$vaccEffEst <- vaccEffEst
    return(parms)
}

testZeros <- function(parmsTmp) {
    tmpCSD <- censSurvDat(parmsTmp)
    casesXgroup <- tmpCSD[,list(cases = sum(infected)), immuneGrp]
    return(0 %in% casesXgroup[,cases])
}

compileStopInfo <- function(minDay, tmp, verbose=0) {
    if(verbose==4) browser()
    out <- data.table(stopDay=minDay
                      , caseCXimmGrpEnd = tmp[immuneGrp==0, sum(infected)]
                      , caseVXimmGrpEnd = tmp[immuneGrp==1, sum(infected)]
                      , hazCXimmGrpEnd = tmp[immuneGrp==0, sum(infected)/sum(perstime)]
                      , hazVXimmGrpEnd = tmp[immuneGrp==1, sum(infected)/sum(perstime)]
                      , ptRatioCVXimmGrpEnd = tmp[immuneGrp==0, sum(perstime)] / tmp[immuneGrp==1, sum(perstime)]
                      )
    out <- as.data.frame(out)
    return(out)
}

getEndResults <- function(parms, bump = T) {
    if(!testZeros(parms)) {
        parmsE <- parms
        parmsE$bump <- F
    }else{
        parmsE <- infBump(parms)
        parmsE$bump <- T
    }
    tmpCSDE <- tmpCSD <- censSurvDat(parms)
    if(testZeros(parms))  tmpCSDE <- censSurvDat(parmsE)
    within(parmsE, {
        if(verbose==2.9) browser()
        vaccEE_ME <- doCoxME(parmsE, tmpCSDE, bump = bump)
        ## vaccEE_GEEclusAR1 <- doGEEclusAR1(clusDat, csd=tmpCSDE, bump = bump)
        ## vaccEE_GLMMclus <- doGLMMclus(parmsE,, csd=tmpCSDE, bump = bump)
        vaccEE_GLMclus <- doGLMclus(parmsE, csd=tmpCSDE, bump = bump)
        vaccEE_GLMFclus <- doGLMFclus(parmsE, csd=tmpCSDE, bump = bump)
        vaccEE_MErelab <- doRelabel(parms, csd=tmpCSD, bump=F, nboot=nboot, verbFreqRelab=10)
        vaccEE_MEboot <- doBoot(parms, csd=tmpCSD, bump=F, nboot=nboot, verbFreqBoot=10)
        vEEs <- list(vaccEE_ME
                     ## , vaccEE_GLMMclus
                     , vaccEE_GLMclus
                     , vaccEE_GLMFclus
                     ## , vaccEE_GEEclusAR1
                     , vaccEE_MEboot
                     , vaccEE_MErelab
                     )
        finMods <- rbindlist(vEEs)
        finInfo <- compileStopInfo(tmp = tmpCSD, minDay=maxInfectDay,  verbose=verbose)
        rm(vaccEE_ME, vaccEE_MEboot, vaccEE_MErelab
           ## , vaccEE_GEEclusAR1
           ## , vaccEE_GLMMclus 
           , vaccEE_GLMFclus , vaccEE_GLMclus
           , vEEs
           )
        return(list(finInfo=finInfo, finMods=finMods))
    })
}

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, returnAll = F,
                       doSeqStops = F, showSeqStops = F, flnm='test', verbose=1, verbFreq=10) {
    set.seed(seed)
    caseXVaccRandGrpList <- caseXPT_ImmuneList <- weeklyAnsList <- list()
    length(caseXVaccRandGrpList) <- length(caseXPT_ImmuneList) <- length(weeklyAnsList) <- N
    finInfo <- finMods <- stopPoints <- data.frame(NULL)
    for(ss in 1:N) {
        if(verbose>0 & (ss %% verbFreq == 0)) print(paste('on',ss,'of',N))
        if(verbose>.5 & (ss %% 1 == 0)) print(paste('on',ss,'of',N))
        if(verbose==2) browser()
        res <- simTrial(parms)
        res <- makeSurvDat(res)
        res <- makeGEEDat(res)
        res <- activeFXN(res)
        res <- getEndResults(res)
        finTmp <- data.frame(sim = ss, res$finMods)
        finMods <- rbind(finMods, finTmp)
        finITmp <- data.frame(sim = ss, res$finInfo)
        finInfo <- rbind(finInfo, finITmp)
        if(doSeqStops) {
            res <- seqStop(res)
            ## if(showSeqStops) {
            ##     resfull <- seqStop(res, fullSeq = T)
            ##     showSeqStop(resfull)
            ## }
            res <- endT(res)
            res <- makeCaseSummary(res)
            stopPt <- as.data.frame(tail(res$weeklyAns,1)) ## active cases by immmune grouping at time of case at end of trial
            stopPt <- with(res, {
                cbind(stopPt
                      , caseCXrandFinA = casesXVaccRandGrp[type=='EVstActive', contCases] ## active cases by vaccination randomization group at final
                      , caseVXrandFinA = casesXVaccRandGrp[type=='EVstActive', vaccCases]
                      , hazCXrandFinA = casesXVaccRandGrp[type=='EVstActive', contCases/contPT]/yearToDays
                      , hazVXrandFinA = casesXVaccRandGrp[type=='EVstActive', vaccCases/vaccPT]/yearToDays
                      , caseCXrandFin = casesXVaccRandGrp[type=='EVst', contCases]         ## total cases by vaccination randomization group at final
                      , caseVXrandFin = casesXVaccRandGrp[type=='EVst', vaccCases]
                      , hazCXrandFin = casesXVaccRandGrp[type=='EVst', contCases/contPT]/yearToDays
                      , hazVXrandFin = casesXVaccRandGrp[type=='EVst', vaccCases/vaccPT]/yearToDays
                      )
            })
            stopPoints <- rbind(stopPoints, stopPt)
        }
        if(returnAll) {
            weeklyAnsList[[ss]] <- as.data.frame(res$weeklyAns)
            caseXVaccRandGrpList[[ss]] <- as.data.frame(res$casesXVaccRandGrp)
            caseXPT_ImmuneList[[ss]] <- as.data.frame(res$casesXPT_Immune)
        }
        rm(res)
        gc()
    }
    if(returnAll)
        return(list(
            stopPoints = stopPoints
            , weeklyAnsList = weeklyAnsList
            , caseXVaccRandGrpList = caseXVaccRandGrpList
            , caseXPT_ImmuneList = caseXPT_ImmuneList
            , finPoint=finPoint
            ))
    if(!returnAll)
        return(list(stopPoints=stopPoints, finMods=finMods, finInfo=finInfo))
}


