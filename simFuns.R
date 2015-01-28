library(survival); library(frailtypack)
# Make a trial population with a given number of clusters of a given size. Put the people in clusters; each cluster has a beta-distributed underlying antibody proportion, which is applied to the people.
makePop <- function(numClus, clusSize){
	pop <- data.table(indiv=as.factor(1:(numClus*clusSize))
		, cluster=as.numeric(gl(n=numClus, k=clusSize))
	)
	return(pop)
}

## Set vaccination time for SWCT by cluster
setSW <- function(p, delayUnit = 7/30, immunoDelay = 21/30) within(p, {
    vaccTime <- delayUnit*(cluster-1)
    immuneTime <- vaccTime + immunoDelay
})

## reparameterize a gamma by mean/var to simulate spatial variation in underlying hazards (change
## later to something more reasonable or based on real data)
reParmRgamma <- function(n, mean, var) {
    theta <- var/mean
    k <- mean/theta
    rgamma(n, shape = k, scale = theta)
}

## Set hazards (monthly units)
setHazs <- function(p, mean, varClus, varIndiv) {
    p[,clusHaz:=NA]
    for(ii in unique(p$cluster)) p <- within(p, { ## set average hazard in cluster
        clusHaz[cluster==ii] <- reParmRgamma(1, mean = mean, var = varClus) 
    })
    p <- within(p, {## individual hazard
        indivHaz <- reParmRgamma(length(indiv), mean = clusHaz, var = varIndiv) 
    })
    return(p)
}

## Simulate SWCT
simSWCT <- function(p, vaccEff = .8, maxInfectTime = 12) {
    vaccRed <- 1 - vaccEff
    p <- within(p, {
        infectTime <- rexp(length(indiv), rate = indivHaz) ## infection pre-vaccination
        notInfectedBeforeVacc <- infectTime > immuneTime ## 
        infectTime[notInfectedBeforeVacc] <- immuneTime[notInfectedBeforeVacc] + rexp(sum(notInfectedBeforeVacc), indivHaz*vaccRed)
        infectTimeTrunc <- infectTime
        infectTimeTrunc[infectTime > maxInfectTime] <- NA
        rm(notInfectedBeforeVacc)
    })
    return(p)
}

## Construct survival data from waiting times
makeSurvDat <- function(p) {
    ## Waiting time up until immune time or vaccination time
    pf <- within(p, {
        startTime <- 0
        endTime <- pmin(immuneTime, infectTime) 
        infected <- as.numeric(infectTime < immuneTime)
        vacc <- 0
    })
    pf <- select(pf, indiv, cluster, vaccTime, immuneTime, startTime, endTime, infected, vacc)
    ## For individuals who experienced vaccination time at risk, tabulate
    pf2 <- within(p[infectTime > immuneTime,], {
        startTime <- immuneTime
        endTime <- infectTime 
        infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
        vacc <- 1
    })
    pf2 <- select(pf2, indiv, cluster, vaccTime, immuneTime, startTime, endTime, infected, vacc)
    pf <- rbind(pf, pf2)
    return(pf)
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

coxPH <- function(csd) { ## take censored survival object and return vacc effectiveness estimates
    ## mod <- coxph(Surv(startTime, endTime, infected) ~ vacc, data=csd) ## without frailty
    modF <- coxph(Surv(startTime, endTime, infected) ~ 
                  vacc + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
                  outer.max=1000, iter.max=10000,
                  data=csd)
    ## 1-summary(mod)$conf.int[,c(1,3:4)], ## without frailty
    vaccEffEst <- 1-summary(modF)$conf.int['vacc',c(1,3:4)] ## gamma frailty
    names(vaccEffEst) <- c('mean','lci','uci')
    return(vaccEffEst)
}

# Make a population and simulate all of the stepped-wedge intervals


# Given a model fit, return a one-tailed P value (goes from 0 for vaccine highly protective to 1 for vaccine highly risky). This is a copy from the other project, which is bad.
oneTailP <- function(m){
	s <- summary(m)
	v <- as.list(s$coef["vacc",])
	if (v$z>0) return(1-v$P/2)
	return(v$P/2)
}

rmodP <- function(people, mform){
	mod <- glmer(mform, data=people,
		family=binomial(link="cloglog")
	)

	if(numTrials<20) print(summary(mod))
 
	return(oneTailP(mod))
}

modP <- function(people, mform){
	mod <- glm(mform, data=people,
		family=binomial(link="cloglog")
	)

	if(numTrials<20) print(summary(mod))
 
	return(oneTailP(mod))
}
