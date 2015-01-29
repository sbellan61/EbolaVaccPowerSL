
## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = 6*30) with(parms, {
    intervalNotStarted <- st[,startDay] > censorDay
    st <- st[!intervalNotStarted,] 
    noInfectionBeforeCensor <- st[,endDay] > censorDay
    st[noInfectionBeforeCensor, infected:=0]
    st[noInfectionBeforeCensor, endDay:=censorDay]
    st[,perstime := (endDay-startDay)]
    st[, active := TRUE] ## for SWCT
    if(trial=='RCT') st[,active :=sum(vacc)>0, by = cluster] ## anyone vaccinated in cluster yet? for RCT analysis
    st <- st[perstime > 0,] 
    return(st)
})

doCoxPH <- function(csd, pkg='coxme', browse=F) { ## take censored survival object and return vacc effectiveness estimates
    if(browse) browser()
    csd <- csd[active==1,]
    if(pkg=='coxph') {
        mod <- try(coxph(Surv(startDay, endDay, infected) ~ vacc + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
                         outer.max=1000, iter.max=10000,
                         data=csd), silent=T)
        vaccEffEst <- 1-summary(mod)$conf.int['vacc',c(1,4:3)] 
        pval <- summary(mod)$coefficients['vacc','p']
        vaccEffEst <- c(vaccEffEst, pval)
    }
    if(pkg=='coxme') {
        mod <- suppressWarnings(coxme(Surv(startDay, endDay, infected) ~ vacc + (1|cluster), data = csd))
        pval <- pnorm(mod$coef/sqrt(vcov(mod)))*2
        vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        vaccEffEst <- c(vaccEffEst, pval)
    }
    ##     mod <- coxph(Surv(startDay, endDay, infected) ~ vacc, data=csd) ## without frailty
    if(inherits(mod, 'try-error')) {print('blah'); browser()}
    names(vaccEffEst) <- c('mean','lci','uci','p')
    return(signif(vaccEffEst,3))
}

doGlmer <- function(csd, bayes=F, browse = F) {## take censored survival object and return vacc effectiveness estimates using bayesian glme
    if(browse) browser()
    csd <- csd[active==1,]
    if(bayes) mod <- bglmer(infected ~ vacc + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    if(!bayes) mod <- glmer(infected ~ vacc + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    vaccRes <- summary(mod)$coefficients['vacc', c('Estimate','Std. Error','Pr(>|z|)')] 
    vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
    names(vaccEffEst) <- c('mean','lci','uci','P')
    return(signif(vaccEffEst,3))
}

summTrial <- function(st) list(summarise(group_by(st, cluster), sum(infected))
                               , summarise(group_by(st, cluster, vacc), sum(infected))
                               , summarise(group_by(st, vacc), sum(infected))
                               )

compileStopInfo <- function(minDay, vaccEffEst, temp) {
    out <- c(stopDay=minDay, vaccEffEst, caseVacc = temp[vacc==1, sum(infected)], caseCont = temp[vacc==0, sum(infected)],
             ptVacc = temp[vacc==1, sum(perstime)], ptCont = temp[vacc==0, sum(perstime)] )
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
    ##     parmsProb <<- parms
    temp <- censSurvDat(parms, midDay)
    vaccEffEst <- doCoxPH(temp) ## converting midDay to days from months
    out <- compileStopInfo(minDay, vaccEffEst, temp)
    if (minDay >= maxDay) return(out)
    pVal <- vaccEffEst['p']
    if(verbose>0) print(signif(vaccEffEst,2)) #paste0('lower 95% of vaccine efficacy at ', signif(midDay,2), ' days =', signif(lciMod,2)))
    goSmaller <- !is.na(pVal) & pVal<.05
    if(goSmaller)
        return(firstStop(parms, minDay, midDay, verbose))
    return(firstStop(parms, midDay+1, maxDay, verbose)) ## output in months
}

seqStop <- function(parms, start = parms$immunoDelay + 14, checkIncrement = 7, verbose = 0, maxDay = 365) {
    trialOngoing <- T
    checkDay <- start
    first <- T
    while(trialOngoing) {
        if(verbose>1) browser()
        temp <- censSurvDat(parms, checkDay)
        vaccEffEst <- try(doCoxPH(temp), silent=T) ## converting midDay to days from months
        if(!inherits(vaccEffEst, 'try-error')) { ## if cox model has enough info to converge check for stopping criteria
            newout <- compileStopInfo(checkDay, vaccEffEst, temp) 
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

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, check=F) {
    set.seed(seed)
    for(ii in 1:N) {
        res <- simTrial(parms)
        stopPoint <- tail(seqStop(res),1)
        if(ii==1) out <- stopPoint else out <- rbind(out, stopPoint)
        if(check) {
            doCoxPH(censSurvDat(res$st, stopPoint$stopDay))
            doCoxPH(censSurvDat(res$st, stopPoint$stopDay+1))
        }
    }
    rownames(out) <- NULL
    return(out)
}

simNwrp <- function(parms=makeParms(), NperCore = 10, check=F, ncores=12) {
    out <- mclapply(1:ncores, simNtrials, N = NperCore, mc.cores = ncores)
    out <- do.call(rbind.data.frame, out)
}
