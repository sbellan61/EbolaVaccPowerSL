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

## Restructure for GEE/GLMM with weekly observations of each cluster.
makeGEEDat <- function(parms, whichDo='pop', verbose=0) within(parms, {
    if(verbose ==1.5) browser()
    popHTmp <- get(paste0(whichDo,'H'))
    popHTmp$immuneDayThink <- popHTmp[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    popHTmp$infectDayRCens <- popHTmp$infectDay
    popHTmp[infectDay==Inf, infectDayRCens := NA]
    popHTmp$immuneGrp <- 0
    popHTmp[day >= immuneDayThink, immuneGrp := 1]
    clusDat <- popHTmp[,list(cases = sum(!is.na(infectDayRCens))), list(cluster, day, immuneGrp)]
    clusDat <- clusDat[order(cluster, day)]
    rm(popHTmp)
})

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

compileStopInfo <- function(minDay, vaccEffEst, tmp, verbose=0) {
    if(verbose==4) browser()
    out <- data.table(stopDay=minDay
                      , mean = vaccEffEst[1,'mean'], lci = vaccEffEst[1,'lci'], uci = vaccEffEst[1,'uci'], p = vaccEffEst[1,'p']
                      , mod = vaccEffEst[1,'mod']
                      , caseCXimmGrpEnd = tmp[immuneGrp==0, sum(infected)]
                      , caseVXimmGrpEnd = tmp[immuneGrp==1, sum(infected)]
                      , hazCXimmGrpEnd = tmp[immuneGrp==0, sum(infected)/sum(perstime)]
                      , hazVXimmGrpEnd = tmp[immuneGrp==1, sum(infected)/sum(perstime)]
                      , ptRatioCVXimmGrpEnd = tmp[immuneGrp==0, sum(perstime)] / tmp[immuneGrp==1, sum(perstime)]
                      )
    out$stopped <- out[, p<.05 & !is.na(p)]
    out$vaccGood <- NA
    out[stopped==T, vaccGood:= lci > 0]
    out <- setcolorder(out, c("stopped", "vaccGood", "stopDay", "mean", "lci", "uci", "p",  'mod'
                              , "caseCXimmGrpEnd", "caseVXimmGrpEnd"
                              , "hazCXimmGrpEnd", "hazVXimmGrpEnd"
                              , "ptRatioCVXimmGrpEnd"
                              ))
    out <- as.data.frame(out)
    return(out)
}

## Check whether stopping point has been reached at intervals
seqStop <- function(parms, start = parms$immunoDelayThink + 14, checkIncrement = 7, minCases = 15,
                    verbose = 0, fullSeq = F, maxDay = parms$maxInfectDay) {
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

getEndResults <- function(parms, verbose = 0) within(parms, {
    tmp <- censSurvDat(parms, maxInfectDay)
    notFitted <- T
    bump <- 0
    tmpElas <- copy(tmp)
    if(verbose==2.9) browser()
    vaccEE_ME <- try(doCoxME(tmpElas, verbose=verbose), silent=T)
    vaccEE_MEboot <- doBoot(tmpElas, nboot=nboot, doFXN = doCoxME, verbose = verbose)
    vaccEE_GEEindivExch <- try(doGEEsurv(tmpElas, verbose = verbose), silent=T)
    vaccEE_GEEclusAR1 <- try(doGEEclusAR1(clusDat, verbose = verbose), silent=T)
    vEEs <- list(vaccEE_ME, vaccEE_MEboot, vaccEE_GEEindivExch, vaccEE_GEEclusAR1)
    stopFin <- rbindlist(lapply(vEEs, compileStopInfo, minDay = maxInfectDay, tmp = tmp, verbose=verbose))
    rm(vaccEE_ME,vaccEE_MEboot, tmp,tmpElas)
})

## while(notFitted) {
## vaccEE_GEEar1 <- try(doGEE(tmpElas, verbose = 3.5, corstr = "ar1"), silent=T)
##    vaccEE_PH <- try(doCoxPH(tmpElas, verbose=verbose), silent=T)
## vaccEE_CL <- try(doGlmer(tmpElas, verbose=verbose), silent=T)
##     if(vaccEE['lci'] > -Inf) notFitted <- F
##     ## pick two individualsin each immune group that weren't infected, and pretend they're infected to do +.5 type analysis in 2x2 table
##     selIndiv <- tmpElas[infectDay==Inf & endDay==168, list(indiv = sample(indiv,1)), immuneGrp]
##     tmpElas[infectDay==Inf & endDay==168 & indiv %in% selIndiv[,indiv], infectDay := endDay] ## perstime still correct
##     tmpElas[endDay==168 & indiv %in% selIndiv[,indiv], infected := 1] 
##     tmpElas[endDay==168 & indiv %in% selIndiv[,indiv], ]
##     bump <- bump+1
## }
## vaccEE <- c(vaccEE, bump=bump)


showSeqStop <- function(resfull, flnm= NULL, ...) {
    with(resfull, {
        if(!is.null(flnm)) pdf(file.path('Figures', paste0(flnm,'.pdf')), ...)
        par(lwd=1.5, mar = c(3,6,1,5), mgp = c(4,1,0), mfrow = c(2,2), oma = c(2,0,0,0))
        plot(0,0, type = 'n', xlim = range(popH$day)/7, ylim = c(0, mu*7), las= 1, bty = 'n',
             ylab = 'weekly hazard', xlab = 'day')
        lines(weeklyAns[,list(stopDay/7, hazCXimmGrpEnd)], col = 'red')
        lines(weeklyAns[,list(stopDay/7, hazVXimmGrpEnd)], col = 'black')
        legend('topleft', c('vacc','cont','total','P value'), col = c('black', 'red','dark green','purple'),  bty = 'n', lwd=2,bg='white')
        ## person-year of observation
        plot(0,0, type = 'n', xlim = range(popH$day)/7, ylim = c(0, yearToDays*max(weeklyAns[, caseCXimmGrpEnd/hazCXimmGrpEnd],na.rm=T)), las= 1, bty = 'n',
             ylab = 'person-years', xlab = 'day')
        lines(weeklyAns[,list(stopDay/7, yearToDays*caseCXimmGrpEnd/hazCXimmGrpEnd)], col = 'red')
        lines(weeklyAns[,list(stopDay/7, yearToDays*caseVXimmGrpEnd/hazVXimmGrpEnd)], col = 'black', lty = 2)
        legend('topleft', c('vacc','cont','P value'), col = c('black', 'red','purple'),  bty = 'n', lwd=2,bg='white')
        ## number of cases
        plot(0,0, type = 'n', xlim = range(popH$day)/7, ylim = c(0, max(weeklyAns[, list(caseVXimmGrpEnd,caseCXimmGrpEnd)])), las= 1, bty = 'n',
             ylab = 'cases', xlab = 'day')
        lines(weeklyAns[,list(stopDay/7, caseCXimmGrpEnd)], col = 'red')
        lines(weeklyAns[,list(stopDay/7, caseVXimmGrpEnd)], col = 'black')
        lines(weeklyAns[,list(stopDay/7, caseVXimmGrpEnd+caseCXimmGrpEnd)], col = 'dark green')
        par(new=T)
        plot(weeklyAns[,list(stopDay, p)], col = 'purple', lty = 1, axes = F, ylab='', xlab='', type='l', ylim = c(0,1))
        axis(4, at = seq(0, 1, by = .05), lab = NA)
        axis(4, at = seq(0, 1, by = .1), las = 1)
        abline(h=.05, lty = 2)
        mtext('p value', 4, 3)
        ## vaccine efficacy estimate
        plot(0,0, type = 'n', xlim = range(popH$day)/7,
             ylim = c(-1, 1), las= 1, bty = 'n', ylab = 'vaccine efficacy', xlab = 'day')
        nmisg <- !weeklyAns[, is.na(lci) | is.na(uci)]
        polygon(c(weeklyAns$stopDay[nmisg], rev(weeklyAns$stopDay[nmisg]))/7,
                c(weeklyAns$lci[nmisg], rev(weeklyAns$uci[nmisg])), col = 'gray', border = NA)
        lines(weeklyAns[,list(stopDay/7, mean)], col = 'black')
        abline(h=0, lty = 2)
        mtext('week of trial', 1, 1, T)
        if(!is.null(flnm)) graphics.off()
    })
}

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, returnAll = F,
                       doSeqStops = F, showSeqStops = F, flnm='test', verbose=1, verbFreq=10) {
    set.seed(seed)
    caseXVaccRandGrpList <- caseXPT_ImmuneList <- weeklyAnsList <- list()
    length(caseXVaccRandGrpList) <- length(caseXPT_ImmuneList) <- length(weeklyAnsList) <- N
    finPoint <- stopPoints <- data.frame(NULL)
    if(showSeqStops) pdf(paste0(flnm, '.pdf'), w = 8, h = 6)
    for(ss in 1:N) {
        if(verbose>0 & (ss %% verbFreq == 0)) print(paste('on',ss,'of',N))
        if(verbose==2) browser()
        ## pseed <- .Random.seed ## for debugging
        ## save(pseed, file = paste0(seed, '-seed.Rdata'))
        res <- simTrial(parms)
        res <- makeSurvDat(res)
        res <- makeGEEDat(res, verbose=verbose)
        res <- activeFXN(res)
        res <- getEndResults(res, verbose=verbose)
        finTmp <- data.frame(sim = ss, res$stopFin)
        finPoint <- rbind(finPoint, finTmp)
        if(doSeqStops) {
            res <- seqStop(res, verbose = 3)
            if(showSeqStops) {
                resfull <- seqStop(res, fullSeq = T)
                showSeqStop(resfull)
            }
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
    if(showSeqStops) graphics.off()
    if(returnAll)
        return(list(
            stopPoints = stopPoints
            , weeklyAnsList = weeklyAnsList
            , caseXVaccRandGrpList = caseXVaccRandGrpList
            , caseXPT_ImmuneList = caseXPT_ImmuneList
            , finPoint=finPoint
            ))
    if(!returnAll)
        return(list(stopPoints=stopPoints, finPoint=finPoint))
}


