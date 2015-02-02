## Construct survival data from waiting times
makeSurvDat <- function(parms) within(parms, {
    pop$immuneDayThink <- pop[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    ## pre-immunity table
    stPre <- copy(pop) # st = survival table
    stPre$startDay <- 0
    stPre$endDay <- stPre[, pmin(immuneDayThink, infectDay)]
    stPre$infected <- stPre[ ,as.numeric(infectDay <= immuneDayThink)]
    stPre$immuneGrp <-  0     ## immuneGrp is variable used for analysis, not omnietient knowledge of vaccination/immune status
    stPre <- stPre[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    ## post-immunity table
    stPost <- copy(pop)[infectDay > immuneDayThink,]
    stPost$startDay <- stPost[,immuneDayThink]
    stPost$endDay   <-  stPost[,infectDay]
    stPost$infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
    stPost$immuneGrp <- 1
    stPost <- stPost[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    st <- rbind(stPre, stPost) ## combine tables
    rm(stPre, stPost)
}) ## careful with modifying parms, st depends on analysis a bit too (immunoDelayThink), so we can have different st for same popH

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = 6*30) with(parms, {
    intervalNotStarted <- st[,startDay] > censorDay
    st <- st[!intervalNotStarted,] 
    noInfectionBeforeCensor <- st[,endDay] > censorDay
    st[noInfectionBeforeCensor, infected:=0]
    st[noInfectionBeforeCensor, endDay:=censorDay]
    st[,perstime := (endDay-startDay)]
    if(trial=='SWCT') ## for SWCT always include all clusters in analysis
        st[, active := TRUE] 
    if(trial %in% c('RCT','FRCT')) ## anyone considered past vaccination refractory in cluster yet? for RCT analysis
        st[,active :=sum(immuneGrp)>0, by = cluster] 
    if(trial=='CRCT') ## for CRCT, only include clusters within a pair that has a cluster that is considered past vaccine refractory
        st[, active := sum(immuneGrp)>0, by = pair] 
    st <- st[perstime > 0,] 
    return(st)
})
## p1 <- simTrial(makeParms(small=F))
## s1 <- makeSurvDat(p1)
## censSurvDat(s1, 72)[,sum(active)]

doCoxPH <- function(csd, pkg='coxme', browse=F) { ## take censored survival object and return vacc effectiveness estimates
    if(browse) browser()
    csd <- csd[active==1,]
    if(pkg=='coxph') {
        mod <- try(coxph(Surv(startDay, endDay, infected) ~ immuneGrp + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
                         outer.max=1000, iter.max=10000,
                         data=csd), silent=T)
        vaccEffEst <- 1-summary(mod)$conf.int['immuneGrp',c(1,4:3)] 
        pval <- summary(mod)$coefficients['immuneGrp','p']
        vaccEffEst <- c(vaccEffEst, pval)
    }
    if(pkg=='coxme') {
        mod <- suppressWarnings(coxme(Surv(startDay, endDay, infected) ~ immuneGrp + (1|cluster), data = csd))
        pval <- pnorm(mod$coef/sqrt(vcov(mod)))*2
        vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        vaccEffEst <- c(vaccEffEst, pval)
    }
    ##     mod <- coxph(Surv(startDay, endDay, infected) ~ immuneGrp, data=csd) ## without frailty
    if(inherits(mod, 'try-error')) {print('blah'); browser()}
    names(vaccEffEst) <- c('mean','lci','uci','p')
    return(signif(vaccEffEst,3))
}

doGlmer <- function(csd, bayes=F, browse = F) {## take censored survival object and return vacc effectiveness estimates using bayesian glme
    if(browse) browser()
    csd <- csd[active==1,]
    if(bayes) mod <- bglmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    if(!bayes) mod <- glmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    vaccRes <- summary(mod)$coefficients['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')] 
    vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
    names(vaccEffEst) <- c('mean','lci','uci','P')
    return(signif(vaccEffEst,3))
}

summTrial <- function(st) list(summarise(group_by(st, cluster), sum(infected))
                               , summarise(group_by(st, cluster, immuneGrp), sum(infected))
                               , summarise(group_by(st, immuneGrp), sum(infected))
                               )

compileStopInfo <- function(minDay, vaccEffEst, tmp) {
    out <- c(stopDay=minDay, vaccEffEst, caseVacc = tmp[immuneGrp==1, sum(infected)], caseCont = tmp[immuneGrp==0, sum(infected)],
             ptVacc = tmp[immuneGrp==1, sum(perstime)], ptCont = tmp[immuneGrp==0, sum(perstime)] )
    out <- c(out, hazVacc = as.numeric(out['caseVacc']/out['ptVacc']/yearToDays),
             hazCont = as.numeric(out['caseCont']/out['ptCont']/yearToDays))
    return(out)
}

## Do a binary search for the number of infections before the stopping point is reached: this is
## assumed to be when 95% CI of vaccine efficacy goes above 0
firstStop <- function(parms, minDay=min(parms$pop$immuneDay) + 30, maxDay=365, verbose = 0) { ## using days to facilitate easier rounding
    if(verbose>=2) browser()
    midDay <- floor((minDay+maxDay)/2) ## floor to days
    ## doCoxPH(censSurvDat(parms$st, midDay), T)
    tmp <- censSurvDat(parms, midDay)
    vaccEffEst <- doCoxPH(tmp) ## converting midDay to days from months
    out <- compileStopInfo(minDay, vaccEffEst, tmp)
    if (minDay >= maxDay) return(out)
    pVal <- vaccEffEst['p']
    if(verbose>0) print(signif(vaccEffEst,2)) #paste0('lower 95% of vaccine efficacy at ', signif(midDay,2), ' days =', signif(lciMod,2)))
    goSmaller <- !is.na(pVal) & pVal<.05
    if(goSmaller)
        return(firstStop(parms, minDay, midDay, verbose))
    return(firstStop(parms, midDay+1, maxDay, verbose)) ## output in months
}

## Check whether stopping point has been reached at intervals
seqStop <- function(parms, start = parms$immunoDelay + 14, checkIncrement = 7, verbose = 0, maxDay = 365) {
    trialOngoing <- T
    checkDay <- start
    first <- T
    while(trialOngoing) {
        if(verbose>1) browser()
        tmp <- censSurvDat(parms, checkDay)
        vaccEffEst <- try(doCoxPH(tmp), silent=T) ## converting midDay to days from months
        if(!inherits(vaccEffEst, 'try-error') & !is.nan(vaccEffEst['p'])) { ## if cox model has enough info to converge check for stopping criteria
            newout <- compileStopInfo(checkDay, vaccEffEst, tmp) 
            if(first) out <- newout else out <- rbind(out, newout)
            first <- F
            if(newout['p'] < .05) trialOngoing <- F
            if(checkDay > maxDay) trialOngoing <- F
        }
        checkDay <- checkDay + 7
    }
    rownames(out) <- NULL
    return(out)
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
