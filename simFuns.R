library(blme); library(survival); library(coxme); library(data.table); library(parallel); library(dplyr); 

yearToDays <- 1/365.25
monthToDays <- 1/30
makeParms <- function(
    trial='RCT'
#  , numClus=20, clusSize=300
  , numClus=4, clusSize=4
  , delayUnit = 7 ## logistically imposed interval in between each new cluster receiving vaccination
  , ord = F ## order clusters' receipt of vaccination
  , mu=.1 * yearToDays ## mean hazard in all participants
  , varClus=mu^2 ##  variance in cluster-level hazards for gamma distribution 
  , sdLogIndiv = 1 ## variance of lognormal distribution of individual RR within a hazard (constant over time, i.e. due to job)
  , vaccEff = .6
  , maxInfectDay = 365 ## end of trial
  , immunoDelay = 21 ## delay from vaccination to immunity
  , weeklyDecay=.9, weeklyDecayVar=.05, hazIntUnit = 7 ## control variation in weekly hazard
){
    if(trial=='FRCT') delayUnit <- delayUnit/2 ## rolling out vaccines as quickly as you would if you were vaccinating whole clusters
    return(as.list(environment()))
}

## create a sequence of weekly hazard trajectories to be used as cluster means
hazTraj <- function(parms=makePop()) within(parms, {
    baseClusHaz <- reParmRgamma(numClus, mean = mu, var = varClus) ## gamma distributed baseline hazards
    decs <- rlnorm(numClus, meanlog = log(weeklyDecay^(1/7)), weeklyDecayVar) ## log-normally distributed decay rates (set var = 0 for constant)
    daySeq <- seq(0,maxInfectDay,by=hazIntUnit)
    hazT <- data.table(day = rep(daySeq, each = numClus), cluster = rep(1:numClus, length(daySeq)), clusHaz = 0)
    cHind <- which(names(hazT)=='clusHaz')
    ## mean cluster hazard trajectory
    for(ii in 1:numClus) hazT[which(hazT[,cluster]==ii), clusHaz := baseClusHaz[ii]*decs[ii]^day]
    ## give every individual a lognormally distributed relative risk
    pop$indivRR <- rlnorm(numClus*clusSize, meanlog = 0, sdlog = sdLogIndiv)
    ## create popH which has weekly hazards for all individuals
    popH <- pop[rep(1:nrow(pop), length(daySeq))]
    popH[, day := rep(daySeq, each=nrow(pop))]
    for(dd in daySeq) for(ii in 1:numClus) popH[day==dd & cluster==ii, clusHaz := hazT[day==dd & cluster==ii, clusHaz]]
    popH[, indivHaz := clusHaz*indivRR]
    rm(ii, cHind, daySeq, decs, baseClusHaz, dd, hazT)
})
## hazTraj(makeParms(weeklyDecay=1, weeklyDecayVar=0))[cluster==1,]
## hazTraj(makeParms(weeklyDecay=.9, weeklyDecayVar=.05))[cluster==1,]

hazTraj()

## Make a trial population with a given number of clusters of a given size. Put the people in
## clusters, give them individual IDs and also id # within cluster
makePop <- function(parms=makeParms()) within(parms, {
    pop <- data.table(indiv=as.factor(1:(numClus*clusSize))
                      , cluster=as.numeric(gl(n=numClus, k=clusSize))
                      , idByClus = rep(1:clusSize, numClus)
                      )
})

## reparameterize a gamma by mean/var to simulate spatial variation in underlying hazards (change
## later to something more reasonable or based on real data)
reParmRgamma <- function(n, mean, var) {
    theta <- var/mean
    k <- mean/theta
    rgamma(n, shape = k, scale = theta)
}

## Set hazards (monthly units)
setHazs <- function(parms) within(parms, {
    pop$indivHaz <- pop$clusHaz<- numeric(nrow(pop)) ## set average hazard in cluster
    cHind <- which(names(pop)=='clusHaz')
    iHind <- which(names(pop)=='indivHaz')
    if(exists('hazT'))
        for(ii in 1:numClus) set(pop, i=which(pop[,cluster]==ii), cHind, reParmRgamma(1, mean = mu, var = varClus))
    else
        for(ii in 1:numClus) set(pop, i=which(pop[,cluster]==ii), cHind, reParmRgamma(1, mean = mu, var = varClus))
    pop[,indivHaz:= reParmRgamma(length(indiv), mean = clusHaz, var = varIndiv)]
    pop[indivHaz<10e-5, indivHaz:=10e-5] ## can't have zero hazard in rexp so make it very small
    rm(cHind, iHind,ii)
})

## Set vaccination time for SWCT by cluster
setSWCTvaccDays <- function(pop, delayUnit = 7, immunoDelay = 21, ord=F) {
    if(ord) { ## if vaccinating highest incidence clusters first (i.e. *NOT* SW randomization of vaccination sequence)
        pop[rev(order(clusHaz))] 
        pp[,indiv:=1:nrow(pp)] ## reset indices so they're ordered again by vaccination sequence
        pp[,cluster:=as.numeric(gl(n=numClus, k=clusSize))]
    }
    pop[, vaccDay := delayUnit*(cluster-1)]
    pop[, immuneDay := vaccDay + immunoDelay]
    return(pop)
}

## Set vaccination time for RCT assuming same speed rollout as SWCT
setRCTvaccDays <- function(pop, delayUnit = 7, immunoDelay = 21, clusSize) {
    if(ord) { ## if vaccinating *HALF* of highest incidence clusters first
        pop[rev(order(clusHaz))] 
        pp[,indiv:=1:nrow(pp)] ## reset indices so they're ordered again by vaccination sequence
        pp[,cluster:=as.numeric(gl(n=numClus, k=clusSize))]
    }
    pop[idByClus <= clusSize/2, vaccDay := Inf]
    pop[idByClus > clusSize/2, vaccDay := delayUnit*(cluster-1)]
    pop[, immuneDay := vaccDay + immunoDelay]
    return(pop)
}

## Set vaccination time for CRCT assuming same speed rollout as SWCT (1 cluster per week)
setCRCTvaccDays <- function(pop, delayUnit = 7, immunoDelay = 21, numClus, clusSize, ord=F) {
    if(ord) { ## if matching clusters on incidence, then randomizing within pairs, then vaccinating highest incidence first
        if(numClus %% 2 == 1) stop("need even # of clusters")
        clusIncRank <- matrix(rev(order(pop[idByClus==1, clusHaz])), nc = 2, byrow=T) ## order clusters by mean hazard, in pairs (each row)
        whichVacc <- sample(1:2, nrow(clusIncRank), replace=T) ## which of each pair to vaccinate
        for(ii in 1:nrow(clusIncRank)) if(whichVacc[ii]==2) clusIncRank[ii,1:2] <- clusIncRank[ii,2:1] ## make first column vaccinated group
        clusIncRank <- c(clusIncRank[,1], clusIncRank[,2]) ## reorder of clusters
        reord <- order(order(rep(clusIncRank, each = clusSize)))
        pop <- pop[reord,] ## reorder
        pop[,indiv:=1:nrow(pop)] ## reset indices so they're ordered again by vaccination sequence
        pop[,cluster:=as.numeric(gl(n=numClus, k=clusSize))]
    }
    pop[cluster <= numClus/2 , vaccDay := delayUnit*(cluster-1)]
    pop[cluster > numClus/2 , vaccDay := Inf]
    pop[, immuneDay := vaccDay + immunoDelay]
    return(pop)
}
## To check ordering works
## simTrial(makeParms('CRCT', numClus=6, clusSize=4, ord=T))$pop[idByClus==1, list(clusHaz*100, rev(order(clusHaz)))]

## Simulate infections
simInfection <- function(pop, vaccEff = .8, maxInfectDay = 12*30) {
    vaccRed <- 1 - vaccEff
    ## infection pre-vaccination
    pop[, infectDay := rexp(length(indiv), rate = indivHaz)] 
    ## infction post-vaccination
    pop[infectDay > immuneDay, infectDay := immuneDay + rexp(length(infectDay), indivHaz*vaccRed)]
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
#    st <- select(st, indiv, cluster, idByClus, vaccDay, immuneDay, startDay, endDay, infected, vacc)
    st <- st[,list(indiv, cluster, idByClus, vaccDay, immuneDay, startDay, endDay, infected, vacc)]
    ## For individuals who experienced vaccination time at risk, tabulate
    st2 <- copy(pop)[infectDay > immuneDay,]
    st2[,startDay:=immuneDay]
    st2[, endDay := infectDay]
    st2[, infected := 1] ## everyone gets infected eventually, but will truncate this in a separate function
    st2[, vacc := 1]
    #st2 <- select(st2, indiv, cluster, idByClus, vaccDay, immuneDay, startDay, endDay, infected, vacc)
    st2 <- st2[,list(indiv, cluster, idByClus, vaccDay, immuneDay, startDay, endDay, infected, vacc)]
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
simTrial <- function(parms=makeParms(), browse = F) suppressWarnings({
    if(browse) browser()
    parms$pop <- do.call(makePop, args=subsArgs(parms, makePop))
    parms <- addDefArgs(parms, makePop)
    parms$pop <- do.call(setHazs, args=subsArgs(parms, setHazs))
    parms <- addDefArgs(parms, setHazs)
    if(parms$trial=='SWCT') {
        parms$pop <- do.call(setSWCTvaccDays, args=subsArgs(parms, setSWCTvaccDays))
        parms <- addDefArgs(parms, setSWCTvaccDays) }
    if(parms$trial %in% c('RCT','FRCT')) {
        parms$pop <- do.call(setRCTvaccDays, args=subsArgs(parms, setRCTvaccDays))
        parms <- addDefArgs(parms, setRCTvaccDays)
    }
    if(parms$trial=='CRCT') {
        parms$pop <- do.call(setCRCTvaccDays, args=subsArgs(parms, setCRCTvaccDays))
        parms <- addDefArgs(parms, setCRCTvaccDays)
    }
    parms$pop <- do.call(simInfection, args=subsArgs(parms, simInfection))
    parms <- addDefArgs(parms, simInfection)
    parms$st <- do.call(makeSurvDat, args=subsArgs(parms, makeSurvDat))
    parms <- addDefArgs(parms, makeSurvDat)
    return(parms)
})
