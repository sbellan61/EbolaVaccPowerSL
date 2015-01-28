library(blmer); library(survival); library(data.table); library(dplyr); library(parallel); 
## Make a trial population with a given number of clusters of a given size. Put the people in
## clusters, give them individual IDs and also id # within cluster
makePop <- function(numClus=20, clusSize=300){
    pop <- data.table(indiv=as.factor(1:(numClus*clusSize))
                      , cluster=as.numeric(gl(n=numClus, k=clusSize))
                      , idByClus = rep(1:clusSize, numClus)
                      )
    return(pop)
}

## Set vaccination time for SWCT by cluster
setSWCTvaccDays <- function(pop, delayUnit = 7, immunoDelay = 21) {
    pop[, vaccDay := delayUnit*(cluster-1)]
    pop[, immuneDay := vaccDay + immunoDelay]
    return(pop)
}

## Set vaccination time for RCT assuming same speed rollout as SWCT
setRCTvaccDays <- function(pop, delayUnit = 7, immunoDelay = 21, clusSize) {
    pop[idByClus <= clusSize/2, vaccDay := Inf]
    pop[idByClus > clusSize/2, vaccDay := delayUnit*(cluster-1)]
    pop[, immuneDay := vaccDay + immunoDelay]
    return(pop)
}

## reparameterize a gamma by mean/var to simulate spatial variation in underlying hazards (change
## later to something more reasonable or based on real data)
reParmRgamma <- function(n, mean, var) {
    theta <- var/mean
    k <- mean/theta
    rgamma(n, shape = k, scale = theta)
}

## Set hazards (monthly units)
setHazs <- function(pop, mu, varClus, varIndiv) {
    pop$indivHaz <- pop$clusHaz<- numeric(nrow(pop)) ## set average hazard in cluster
    cHind <- which(names(pop)=='clusHaz')
    iHind <- which(names(pop)=='indivHaz')
    for(ii in unique(pop$cluster)) set(pop, i=which(pop[,cluster]==ii), cHind, reParmRgamma(1, mean = mu, var = varClus))
    pop[,indivHaz:= reParmRgamma(length(indiv), mean = clusHaz, var = varIndiv)]
    pop[indivHaz<10e-5, indivHaz:=10e-5] ## can't have zero hazard in rexp so make it very small
    return(pop)
}

## Simulate infections
simInfection <- function(pop, vaccEff = .8, maxInfectDay = 12*30) {
    vaccRed <- 1 - vaccEff
    ## infection pre-vaccination
    pop[, infectDay := rexp(length(indiv), rate = indivHaz)] 
    notInfectedBeforeVacc <- with(pop, infectDay > immuneDay) 
    ## infction post-vaccination
    if(sum(notInfectedBeforeVacc)>0) pop[notInfectedBeforeVacc, infectDay := immuneDay + rexp(sum(notInfectedBeforeVacc), indivHaz*vaccRed)]
    pop[, infectDayTrunc := infectDay]
    pop[infectDay > maxInfectDay, infectDayTrunc := NA]
    return(pop)
}

## Construct survival data from waiting times
makeSurvDat <- function(pop) {
    ## Waiting time up until immune time or vaccination time
    st <- copy(pop)
    st[,startDay:=0]
    st[,endDay:=pmin(immuneDay, infectDay)]
    st[,infected:=as.numeric(infectDay < immuneDay)]
    st[,vacc:=0]
    st <- select(st, indiv, cluster, idByClus, vaccDay, immuneDay, startDay, endDay, infected, vacc)
    ## For individuals who experienced vaccination time at risk, tabulate
    st2 <- copy(pop)[infectDay > immuneDay,]
    st2[,startDay:=immuneDay]
    st2[, endDay := infectDay]
    st2[, infected := 1] ## everyone gets infected eventually, but will truncate this in a separate function
    st2[, vacc := 1]
    st2 <- select(st2, indiv, cluster, idByClus, vaccDay, immuneDay, startDay, endDay, infected, vacc)
    st <- rbind(st, st2)
    return(st)
}

subsArgs <- function(parms, fxn) parms[names(parms) %in% names(formals(fxn))] ## get parameters necessary for a fxn
addDefArgs <- function(parms, fxn) { ## add default arguments to parameter list if not specified in parameter list
    pnms <- names(parms)
    fnms <- names(formals(fxn))
    for(nn in fnms[!fnms %in% pnms]) parms[nn] <- formals(fxn)[nn]
    parms[!names(parms) %in% names(formals(fxn))]
    return(parms)
}

## simulate whole trial and return with all parameters used
simTrial <- function(parms=makeParms(), browse = F) {
    if(browse) browser()
    parms$pop <- do.call(makePop, args=subsArgs(parms, makePop))
    parms <- addDefArgs(parms, makePop)
    if(parms$trial=='SWCT') {
        parms$pop <- do.call(setSWCTvaccDays, args=subsArgs(parms, setSWCTvaccDays))
        parms <- addDefArgs(parms, setSWCTvaccDays) }
    if(parms$trial=='RCT') {
        parms$pop <- do.call(setRCTvaccDays, args=subsArgs(parms, setRCTvaccDays))
        parms <- addDefArgs(parms, setRCTvaccDays)
    }
    parms$pop <- do.call(setHazs, args=subsArgs(parms, setHazs))
    parms <- addDefArgs(parms, setHazs)
    parms$pop <- do.call(simInfection, args=subsArgs(parms, simInfection))
    parms <- addDefArgs(parms, simInfection)
    parms$st <- do.call(makeSurvDat, args=subsArgs(parms, makeSurvDat))
    parms <- addDefArgs(parms, makeSurvDat)
    return(parms)
}

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = 6*30) with(parms, {
    intervalNotStarted <- st[,startDay] > censorDay
    st <- st[!intervalNotStarted,] 
    noInfectionBeforeCensor <- st[,endDay] > censorDay
    st[noInfectionBeforeCensor, infected:=0]
    st[noInfectionBeforeCensor, endDay:=censorDay]
    st[,perstime := (endDay-startDay)]
    st[,active :=sum(vacc)>0, by = cluster] ## anyone vaccinated in cluster yet? for RCT analysis
    st <- st[perstime > 0,] 
    return(st)
})

doCoxPH <- function(csd, frail=T,browse=F) { ## take censored survival object and return vacc effectiveness estimates
    ##csd <- csd[clusActive==1,]
    if(browse) browser()
    if(frail) {    mod <- coxph(Surv(startDay, endDay, infected) ~ 
                                vacc + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
                                outer.max=1000, iter.max=10000,
                                data=csd)
               }else{
                   mod <- coxph(Surv(startDay, endDay, infected) ~ vacc, data=csd) ## without frailty
               }
    vaccEffEst <- 1-summary(mod)$conf.int['vacc',c(1,4:3)] ## gamma frailty
    names(vaccEffEst) <- c('mean','lci','uci')
    pval <- ifelse(frail, summary(mod)$coefficients['vacc','p'], summary(mod)$coefficients['Pr(>|z|)'])
    vaccEffEst <- c(vaccEffEst, P = pval)
                                        #    if(is.na(vaccEffEst['P'])) vaccEffEst <- doCoxPH(csd, frail=F)
    return(signif(vaccEffEst,3))
}

doGlmer <- function(csd, bayes=F, browse = F) {## take censored survival object and return vacc effectiveness estimates using bayesian glme
    if(browse) browser()
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

## Do a binary search for the number of infections before the stopping point is reached: this is
## assumed to be when 95% CI of vaccine efficacy goes above 0
firstStop <- function(parms, minDay=min(parms$pop$immuneDay) + 30, maxDay=365, verbose = 0) { ## using days to facilitate easier rounding
if(verbose>=2) browser()
    midDay <- floor((minDay+maxDay)/2) ## floor to days
## doCoxPH(censSurvDat(parms$st, midDay), T)
##     parmsProb <<- parms
    vaccEffEst <- doCoxPH(censSurvDat(parms$st, midDay)) ## converting midDay to days from months
    if (minDay >= maxDay) return(c(stopDay=minDay, vaccEffEst))
    pVal <- vaccEffEst['P']
    if(verbose>0) print(signif(vaccEffEst,2)) #paste0('lower 95% of vaccine efficacy at ', signif(midDay,2), ' days =', signif(lciMod,2)))
    if(pVal<0.05)
        return(firstStop(parms, minDay, midDay, verbose))
    return(firstStop(parms, midDay+1, maxDay, verbose)) ## output in months
}

casesInTrial <- function(parms, maxDayCaseDay = 6*30) sum(with(parms$pop, infectDay < maxDayCaseDay))


yearToDays <- 1/365.25
monthToDays <- 1/30
makeParms <- function(trial='RCT', 
                      mu=.1 / 365.25, varClus=mu^2/2, varIndiv = mu^2/8,  ## hazards
                      vaccEff = .6, maxInfectDay = 12*30,
                      numClus=20, clusSize=300){
    list(mu=mu, varClus=varClus, varIndiv = varIndiv, trial=trial, vaccEff = vaccEff, maxInfectDay=maxInfectDay,
         numClus=numClus, clusSize=clusSize)
}

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, check=F) {
    set.seed(seed)
    for(ii in 1:N) {
        res <- simTrial(parms)
        stopPoint <- firstStop(res)
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
