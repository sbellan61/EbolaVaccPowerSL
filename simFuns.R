library(blme); library(survival); library(coxme); library(data.table); library(parallel); library(dplyr); 

yearToDays <- 1/365.25
monthToDays <- 1/30
makeParms <- function(
    trial='RCT'
    , numClus=20, clusSize=300
    , delayUnit = 7 ## logistically imposed interval in between each new cluster receiving vaccination
    , ord = 'none' ## order clusters' receipt of vaccination ('none', by baseline visit 'BL', by time-updated 'TU' last interval incidence)
    , mu=.1 * yearToDays ## mean hazard in all participants
    , varClus=mu^2 ##  variance in cluster-level hazards for gamma distribution 
    , sdLogIndiv = 1 ## variance of lognormal distribution of individual RR within a hazard (constant over time, i.e. due to job)
    , vaccEff = .6
    , maxInfectDay = 35 ## end of trial
    , immunoDelay = 21 ## delay from vaccination to immunity
    , immunoDelayThink = immunoDelay ## delay from vaccination to immunity used in analysis (realistically would be unknown)
    , weeklyDecay=.9, weeklyDecayVar=.05 ## log-normally distributed incidence decay rates (set var = 0 for constant)
    , hazIntUnit = 7 ## interval between discrete changes in hazard
    , reordLag = 0 ## how long ago's hazard to use when deciding this week's time-updated vaccination sequence
    , small=F ## do a small trial for illustration
    ){
    if(small) {
        numClus <- 4
        clusSize <- 3
    }
    if(trial=='FRCT') delayUnit <- delayUnit/2 ## rolling out vaccines as quickly as you would if you were vaccinating whole clusters
    return(as.list(environment()))
}

## Make a trial population with a given number of clusters of a given size. Put the people in
## clusters, give them individual IDs and also id # within cluster
makePop <- function(parms=makeParms()) within(parms, {
    pop <- data.table(indiv=as.factor(1:(numClus*clusSize))
                      , cluster=as.numeric(gl(n=numClus, k=clusSize))
                      , idByClus = rep(1:clusSize, numClus)
                      )
})

## Set cluster- and individual-level hazards, with cluster means changing over time and individual
## RR around cluster mean constant
setHazs <- function(parms=makePop()) within(parms, {
    baseClusHaz <- reParmRgamma(numClus, mean = mu, var = varClus) ## gamma distributed baseline hazards
    decs <- rlnorm(numClus, meanlog = log(weeklyDecay^(1/7)), weeklyDecayVar) 
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
    rm(ii, cHind, decs, baseClusHaz, dd, hazT)
})
## setHazs(makePop(makeParms(weeklyDecay=1, weeklyDecayVar=0)))$popH[cluster==1,]
## setHazs(makePop(makeParms(weeklyDecay=.9, weeklyDecayVar=0)))$popH[cluster==1,]

## reparameterize a gamma by mean/var to simulate spatial variation in underlying hazards (change
## later to something more reasonable or based on real data)
reParmRgamma <- function(n, mean, var) {
    theta <- var/mean
    k <- mean/theta
    rgamma(n, shape = k, scale = theta)
}

## Set vaccination time for SWCT by cluster
setSWCTvaccDays <- function(parms) {
    parms <- reordSWCTorRCT(parms)
    within(parms, {
        popH$vaccDay <- delayUnit*(popH[,cluster]-1)
        popH$immuneDay <- popH[,vaccDay] + immunoDelay
        popH$immuneDayThink <- popH[,vaccDay] + immunoDelayThink ## refractory period for vaccine ASSUMED in analysis
        popH$vacc <- popH[, day>=vaccDay]
        popH$immune <- popH[, day>=immuneDay]
        ## reset pop to refrence data table
        pop <- select(popH[day==0,], indiv, cluster, idByClus, indivRR, vaccDay, immuneDay, immuneDayThink) 
    })
}
p1 <- setHazs(makePop(makeParms(clusSize=2, weeklyDecay=.9, weeklyDecayVar=.2, ord='TU', trial='SWCT',small=T)))
p1 <- setVaccDays(p1)
p1$popH[idByClus==1,list(cluster,clusHaz, day,vacc,immune)]

## Set vaccination time for RCT assuming same speed rollout as SWCT
setRCTvaccDays <- function(parms) within(parms, {
    parms <- reordSWCTorRCT(parms)
    popH$vaccDay <- Inf ## unvaccinated
    popH[idByClus > clusSize/2 , vaccDay := delayUnit*(cluster-1)] ## half get vaccinated
    popH$immuneDay <- popH[,vaccDay] + immunoDelay
    popH$immuneDayThink <- popH[,vaccDay] + immunoDelayThink ## refractory period for vaccine ASSUMED in analysis
    popH$vacc <- popH[, day>=vaccDay]
    popH$immune <- popH[, day>=immuneDay]
    pop <- select(popH[day==0,], indiv, cluster, idByClus, indivRR, vaccDay, immuneDay, immuneDayThink) ## reset pop to refrence data table
})

reordSWCTorRCT <- function(parms) within(parms, {
    if(ord=='BL') { ## if vaccinating highest incidence clusters first (i.e. *NOT* SW randomization of vaccination sequence)
        clusIncRank <- popH[idByClus==1 & day == 0,order(rev(order(clusHaz)))]
        popH[, cluster:=clusIncRank[popH[, cluster]]]
        popH <- arrange(popH, day, cluster)
        popH[,indiv:= rep(1:(numClus*clusSize), length(daySeq))] ## reset indices so they're ordered again by vaccination sequence
        rm(clusIncRank)
    }
    if(ord=='TU') { ## time-updated sequencing, each week vaccinate highest incidence cluster from reordLag days ago
        clusIncRank <- NULL
        for(ii in 1:numClus) { ## for each vaccination day
            dd <- daySeq[ii]
            updatingOrder <- 1:numClus > ii-1 ## i.e. on 3rd day of vaccination, only updating the 3rd vaccination sequence
            currentRank <- popH[idByClus==1 & day == dd - reordLag, rev(order(clusHaz))] ## current cluster hazard ordering
            currentRank <- currentRank[!currentRank %in% clusIncRank] ## remove clusters already vaccinated
            clusIncRank <- c(clusIncRank, currentRank[1])
        }
        clusIncRank <- order(clusIncRank)
        popH[, cluster:=clusIncRank[popH[, cluster]]]
        popH <- arrange(popH, day, cluster)
        popH[,indiv:= rep(1:(numClus*clusSize), length(daySeq))] ## reset indices so they're ordered again by vaccination sequence
    }
})


## Set vaccination time for CRCT assuming same speed rollout as SWCT (1 cluster per week)
setCRCTvaccDays <- function(parms) within(parms, {
    if(ord=='BL') { ## if matching clusters on incidence, then randomizing within pairs, then vaccinating highest incidence first
        if(numClus %% 2 == 1) stop("need even # of clusters")
        ## order clusters by mean hazard, in pairs (each row)
        clusIncRank <- matrix(rev(order(popH[day==0 & idByClus==1, clusHaz])), nc = 2, byrow=T) 
        whichVacc <- sample(1:2, nrow(clusIncRank), replace=T) ## which of each pair to vaccinate
        for(ii in 1:nrow(clusIncRank)) if(whichVacc[ii]==2) clusIncRank[ii,1:2] <- clusIncRank[ii,2:1] ## make first column vaccinated group
        clusIncRank <- c(clusIncRank[,1], clusIncRank[,2]) ## reorder of clusters
        clusIncRank <- order(clusIncRank)
        popH[, cluster:=clusIncRank[popH[, cluster]]]
        popH <- arrange(popH, day, cluster)
        popH[,indiv:= rep(1:(numClus*clusSize), length(daySeq))] ## reset indices so they're ordered again by vaccination sequence
        rm(clusIncRank)
    }
    popH$vaccDay <- Inf
    popH[cluster <= numClus/2, vaccDay := delayUnit*(cluster-1)]
    popH$immuneDay <- popH[,vaccDay] + immunoDelay ## true refractory period for vaccine
    popH$immuneDayThink <- popH[,vaccDay] + immunoDelayThink ## refractory period for vaccine ASSUMED in analysis
    popH$vacc <- popH[, day>=vaccDay]
    popH$immune <- popH[, day>=immuneDay]
    pop <- select(popH[day==0,], indiv, cluster, idByClus, indivRR, vaccDay, immuneDay, immuneDayThink) ## reset pop to refrence data table
})
## To check ordering works
## p1 <- setHazs(makePop(makeParms(clusSize=2, weeklyDecay=.9, weeklyDecayVar=0, ord='BL')))
## setVaccDays(p1)$popH[,list(cluster,clusHaz, day,vacc,immune)]

setVaccDays <- function(parms) { ## wrapper around other functions above
    setVaccFXN <- get(paste0('set',parms$trial,'vaccDays'))
    setVaccFXN(parms)
}

## Simulate infections
simInfection <- function(parms) within(parms, {
    popH$infectDay <- Inf
    for(dd in daySeq) { ## infection day is beginning of each hazard interval + exponential waiting time
        popH[day==dd, infectDay := dd + rexp(length(indiv), rate = indivHaz*ifelse(immune, 1-vaccEff, 1))] 
        popH[day==dd & infectDay > dd + hazIntUnit, infectDay := Inf] ## reset if it goes into next haard interval
    }
    pop$infectDay <- Inf ## copy infection days to pop, to use in analysis
    pop[indiv == popH[infectDay!=Inf, indiv], infectDay:= popH[infectDay!=Inf, infectDay]]
    rm(dd)
})

## p1 <- setHazs(makePop(makeParms(clusSize=300, numClus=20, weeklyDecay=.9, weeklyDecayVar=0, ord='BL')))
## p1 <- setSWCTvaccDays(p1)
## p1 <- simInfection(p1)
## head(p1$pop[infectDay!=Inf, list(cluster, immuneDay, infectDay)],100)

## Construct survival data from waiting times
makeSurvDat <- function(parms) within(parms, {
    ## pre-immunity table
    stPre <- copy(pop) # st = survival table
    stPre$startDay <- 0
    stPre$endDay <- stPre[, pmin(immuneDayThink, infectDay)]
    stPre$infected <- stPre[ ,as.numeric(infectDay < immuneDayThink)]
    stPre$immuneGrp <-  0     ## immuneGrp is variable used for analysis, not omnietient knowledge of vaccination/immune status
    stPre <- stPre[,list(indiv, cluster, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    ## post-immunity table
    stPost <- copy(pop)[infectDay > immuneDayThink,]
    stPost$startDay <- stPost[,immuneDayThink]
    stPost$endDay   <-  stPost[,infectDay]
    stPost$infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
    stPost$immuneGrp <- 1
    stPost <- stPost[,list(indiv, cluster, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    st <- rbind(stPre, stPost) ## combine tables
    rm(stPre, stPost)
})

## simulate whole trial and return with all parameters used
simTrial <- function(parms=makeParms(), browse = F) {
    if(browse) browser()
    parms <- makePop(parms)
    parms <- setHazs(parms)
    parms <- setVaccDays(parms)
    parms <- simInfection(parms)
    parms <- makeSurvDat(parms)
    return(parms)
}
simTrial()
