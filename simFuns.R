library(blme); library(survival); library(coxme); library(data.table); library(parallel); 
#library(dplyr); 

yearToDays <- 1/365.25
monthToDays <- 1/30
makeParms <- function(trial='RCT', delayUnit = 7,
                      mu=.1 / 365.25, varClus=mu^2, varIndiv = mu^2/8,  ## hazards
                      vaccEff = .6, maxInfectDay = 12*30, immunoDelay = 21,
                      numClus=20, clusSize=300){
    list(mu=mu, varClus=varClus, varIndiv = varIndiv, trial=trial, vaccEff = vaccEff, maxInfectDay=maxInfectDay,
         numClus=numClus, clusSize=clusSize, delayUnit=delayUnit)
}

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

## Set vaccination time for CRCT assuming same speed rollout as SWCT (1 cluster per week)
setCRCTvaccDays <- function(pop, delayUnit = 7, immunoDelay = 21, numClus) {
    pop[cluster <= numClus/2 , vaccDay := delayUnit*(cluster-1)]
    pop[cluster > numClus/2 , vaccDay := Inf]
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

## runFxn <- function(fxn, parms) {
##     do.call(fxn, args=subsArgs(parms, fxn))
##     parms <- addDefArgs(parms, fxn)
##     return(parms)
## }

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
    if(parms$trial=='CRCT') {
        parms$pop <- do.call(setCRCTvaccDays, args=subsArgs(parms, setCRCTvaccDays))
        parms <- addDefArgs(parms, setCRCTvaccDays)
    }
    parms$pop <- do.call(setHazs, args=subsArgs(parms, setHazs))
    parms <- addDefArgs(parms, setHazs)
    parms$pop <- do.call(simInfection, args=subsArgs(parms, simInfection))
    parms <- addDefArgs(parms, simInfection)
    parms$st <- do.call(makeSurvDat, args=subsArgs(parms, makeSurvDat))
    parms <- addDefArgs(parms, makeSurvDat)
    return(parms)
}
