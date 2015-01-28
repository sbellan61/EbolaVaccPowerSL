library(survival); library(frailtypack)
# Make a trial population with a given number of clusters of a given size. Put the people in clusters; each cluster has a beta-distributed underlying antibody proportion, which is applied to the people.
makePop <- function(numClus=20, clusSize=300){
	pop <- data.table(indiv=as.factor(1:(numClus*clusSize))
		, cluster=as.numeric(gl(n=numClus, k=clusSize))
	)
	return(pop)
}

## Set vaccination time for SWCT by cluster
setSW <- function(pop, delayUnit = 7/30, immunoDelay = 21/30) {
    pop[, vaccTime := delayUnit*(cluster-1)]
    pop[, immuneTime := vaccTime + immunoDelay]
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
    return(pop)
}

## Simulate SWCT
simSWCTinfection <- function(pop, vaccEff = .8, maxInfectTime = 12) {
    vaccRed <- 1 - vaccEff
    ## infection pre-vaccination
    pop[, infectTime := rexp(length(indiv), rate = indivHaz)] 
    notInfectedBeforeVacc <- with(pop, infectTime > immuneTime) 
    ## infction post-vaccination
    pop[notInfectedBeforeVacc, infectTime := immuneTime + rexp(sum(notInfectedBeforeVacc), indivHaz*vaccRed)]
    pop[, infectTimeTrunc := infectTime]
    pop[infectTime > maxInfectTime, infectTimeTrunc := NA]
    return(pop)
}

## Construct survival data from waiting times
makeSurvDat <- function(pop) {
    ## Waiting time up until immune time or vaccination time
    st <- copy(pop)
    st[,startTime:=0]
    st[,endTime:=pmin(immuneTime, infectTime)]
    st[,infected:=as.numeric(infectTime < immuneTime)]
    st[,vacc:=0]
    st <- select(st, indiv, cluster, vaccTime, immuneTime, startTime, endTime, infected, vacc)
    st
    ## For individuals who experienced vaccination time at risk, tabulate
    st2 <- copy(pop)[infectTime > immuneTime,]
    st2[,startTime:=immuneTime]
    st2[, endTime := infectTime]
    st2[, infected := 1] ## everyone gets infected eventually, but will truncate this in a separate function
    st2[, vacc := 1]
    st2 <- select(st2, indiv, cluster, vaccTime, immuneTime, startTime, endTime, infected, vacc)
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
## simulate whole SWCT and return with all parameters used
simSWCT <- function(parms=list(mu=.0065, varClus=2e-05, varIndiv = 5e-06)) {
    parms$pop <- do.call(makePop, args=subsArgs(parms, makePop))
    parms <- addDefArgs(parms, makePop)
    parms$pop <- do.call(setSW, args=subsArgs(parms, setSW))
    parms <- addDefArgs(parms, setSW)
    parms$pop <- do.call(setHazs, args=subsArgs(parms, setHazs))
    parms <- addDefArgs(parms, setHazs)
    parms$pop <- do.call(simSWCTinfection, args=subsArgs(parms, simSWCTinfection))
    parms <- addDefArgs(parms, simSWCTinfection)
    parms$st <- do.call(makeSurvDat, args=subsArgs(parms, makeSurvDat))
    parms <- addDefArgs(parms, makeSurvDat)
    return(parms)
}

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(st, censorTime = 6) {
    intervalNotStarted <- st[,startTime] > censorTime
    st <- st[!intervalNotStarted,] 
    noInfectionBeforeCensor <- st[,endTime] > censorTime
    st[noInfectionBeforeCensor, infected:=0]
    st[noInfectionBeforeCensor, endTime:=censorTime]
    return(st)
}

doCoxPH <- function(csd, browse=F) { ## take censored survival object and return vacc effectiveness estimates
    ## mod <- coxph(Surv(startTime, endTime, infected) ~ vacc, data=csd) ## without frailty
    if(browse) browser()
    modF <- coxph(Surv(startTime, endTime, infected) ~ 
                  vacc + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
                  outer.max=1000, iter.max=10000,
                  data=csd)
    ## 1-summary(mod)$conf.int[,c(1,4:3)], ## without frailty
    vaccEffEst <- 1-summary(modF)$conf.int['vacc',c(1,4:3)] ## gamma frailty
    names(vaccEffEst) <- c('mean','lci','uci')
    return(vaccEffEst)
}

## Do a binary search for the number of infections before the stopping point is reached: this is
## assumed to be when 95% CI of vaccine efficacy goes above 0
firstStop <- function(parms, min=7, max=12*30, verbose = 0) { ## using days to facilitate easier rounding
    if (min >= max) return(min)
    mid <- floor((min+max)/2) ## floor to days
    lciMod <- doCoxPH(censSurvDat(parms$st, mid/30))['lci'] ## converting mid to days from months
    if(verbose>0) print(paste0('lower 95% of vaccine efficacy at ', signif(mid/30,2), ' months =', signif(lciMod,2)))
    if(lciMod>=0)
        return(firstStop(parms, min, mid, verbose))
    return(firstStop(parms, mid+1, max, verbose)/30) ## output in months
}

## # Given a model fit, return a one-tailed P value (goes from 0 for vaccine highly protective to 1 for vaccine highly risky). This is a copy from the other project, which is bad.
## oneTailP <- function(m){
## 	s <- summary(m)
## 	v <- as.list(s$coef["vacc",])
## 	if (v$z>0) return(1-v$P/2)
## 	return(v$P/2)
## }

## rmodP <- function(people, mform){
## 	mod <- glmer(mform, data=people,
## 		family=binomial(link="cloglog")
## 	)

## 	if(numTrials<20) print(summary(mod))
 
## 	return(oneTailP(mod))
## }

## modP <- function(people, mform){
## 	mod <- glm(mform, data=people,
## 		family=binomial(link="cloglog")
## 	)

## 	if(numTrials<20) print(summary(mod))
 
## 	return(oneTailP(mod))
## }
