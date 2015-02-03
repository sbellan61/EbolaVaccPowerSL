## Construct survival data from waiting times
makeSurvDat <- function(parms,  whichDo='pop') within(parms, {
    popTmp <- get(whichDo)
    popTmp$immuneDayThink <- popTmp[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    ## pre-immunity table
    stPre <- copy(popTmp) # st = survival table
    stPre$startDay <- 0
    stPre$endDay <- stPre[, pmin(immuneDayThink, infectDay)]
    stPre$infected <- stPre[ ,as.numeric(infectDay <= immuneDayThink)]
    stPre$immuneGrp <-  0     ## immuneGrp is variable used for analysis, not omnietient knowledge of vaccination/immune status
    stPre <- stPre[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    ## post-immunity table
    stPost <- copy(popTmp)[infectDay > immuneDayThink,]
    stPost$startDay <- stPost[,immuneDayThink]
    stPost$endDay   <-  stPost[,infectDay]
    stPost$infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
    stPost$immuneGrp <- 1
    stPost <- stPost[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    nmSt <- paste0(sub('pop','',whichDo),'st') ## makes EVpop into EVst for example
    assign(nmSt, rbind(stPre, stPost)) ## combine tables
    rm(stPre, stPost, popTmp, nmSt)
}) ## careful with modifying parms, st depends on analysis a bit too (immunoDelayThink), so we can have different st for same popH

## Select subset of survival table to analyze
activeFXN <- function(parms, whichDo='st') within(parms, { 
    ## for SWCT or unmatched CRCT always include all clusters in analysis because unvaccinated
    ## clusters are still considered to provide useful information from baseline
    stA <- copy(st)
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
censSurvDat <- function(parms, censorDay = 6*30, whichDo = 'stActive') with(parms, {
    stTmp <- get(whichDo)
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
    out <- c(stopDay=minDay, vaccEffEst, caseVacc = tmp[immuneGrp==1, sum(infected)], caseCont = tmp[immuneGrp==0, sum(infected)],
             ptVacc = tmp[immuneGrp==1, sum(perstime)], ptCont = tmp[immuneGrp==0, sum(perstime)] )
    out <- c(out, hazVacc = as.numeric(out['caseVacc']/out['ptVacc']/yearToDays),
             hazCont = as.numeric(out['caseCont']/out['ptCont']/yearToDays))
    out <- c(out, ptRatio = as.numeric(out['ptCont'] / out['ptVacc']))
    return(out)
}

## ## Do a binary search for the number of infections before the stopping point is reached: this is
## ## assumed to be when 95% CI of vaccine efficacy goes above 0
## firstStop <- function(parms, minDay=min(parms$pop$immuneDayThink) + 30, maxDay=365, verbose = 0) { ## using days to facilitate easier rounding
##     if(verbose>=2) browser()
##     midDay <- floor((minDay+maxDay)/2) ## floor to days
##     tmp <- censSurvDat(parms, midDay)
##     vaccEffEst <- doCoxPH(tmp) ## converting midDay to days from months
##     out <- compileStopInfo(minDay, vaccEffEst, tmp)
##     if (minDay >= maxDay) return(out)
##     pVal <- vaccEffEst['p']
##     if(verbose>0) print(signif(vaccEffEst,2))
##     goSmaller <- !is.na(pVal) & pVal<.05
##     if(goSmaller)
##         return(firstStop(parms, minDay, midDay, verbose))
##     return(firstStop(parms, midDay+1, maxDay, verbose)) ## output in months
## }

## Check whether stopping point has been reached at intervals
seqStop <- function(parms, start = parms$immunoDelayThink + 14, checkIncrement = 7, verbose = 0, maxDay = parms$maxInfectDay) {
    trialOngoing <- T
    checkDay <- start
    first <- T
    while(trialOngoing) {
        if(verbose>1) browser()
        if(verbose>0) print(checkDay)
        tmp <- censSurvDat(parms, checkDay)
        vaccEffEst <- try(doCoxPH(tmp), silent=F) ## converting midDay to days from months
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



casesInTrial <- function(parms, maxDayCaseDay = 6*30) sum(with(parms$pop, infectDay < maxDayCaseDay))

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, check=F, verbose=0) {
    set.seed(seed)
    for(ii in 1:N) {
        if(verbose>1) browser()
        res <- simTrial(parms)
        stopPoint <- tail(seqStop(res),1)
        if(ii==1) out <- stopPoint else out <- rbind(out, stopPoint)
        if(check) {
            doCoxPH(censSurvDat(res$st, stopPoint$stopDay))
            doCoxPH(censSurvDat(res$st, stopPoint$stopDay+1))
        }
    }
    rownames(out) <- NULL
    return(out)#as.data.table(out))
}

simNwrp <- function(parms=makeParms(), NperCore = 10, check=F, ncores=12) {
    out <- mclapply(1:ncores, simNtrials, N = NperCore, mc.cores = ncores)
    out <- do.call(rbind.data.frame, out)
}
