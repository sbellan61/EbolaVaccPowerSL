library(blme); library(survival); library(coxme); library(data.table); library(parallel); library(dplyr); 

yearToDays <- 1/365.25
monthToDays <- 1/30
trialTypes <- c('RCT','FRCT','SWCT','CRCT')
makeParms <- function(
    trial='RCT'
    , numClus=20, clusSize=300
    , delayUnit = 7 ## logistically imposed interval in between each new cluster receiving vaccination
    , ord = 'none' ## order clusters' receipt of vaccination ('none', by baseline visit 'BL', by time-updated 'TU' last interval incidence)
    , mu=.03 * yearToDays ## mean hazard in all participants at baseline
    , varClus=mu^2/2 ##  variance in cluster-level hazards for gamma distribution 
    , sdLogIndiv = 1 ## variance of lognormal distribution of individual RR within a hazard (constant over time, i.e. due to job)
    , vaccEff = .8
    , maxInfectDay = 7*36 ## end of trial (36 weeks default; 9 months)
    , immunoDelay = 21 ## delay from vaccination to immunity
    , immunoDelayThink = immunoDelay ## delay from vaccination to immunity used in analysis (realistically would be unknown)
    , weeklyDecay=.9, weeklyDecayVar=.007 ## log-normally distributed incidence decay rates (set var = 0 for constant)
    , hazIntUnit = 7 ## interval between discrete changes in hazard
    , reordLag = 14 ## how long ago's hazard to use when deciding this week's time-updated vaccination sequence
    , includeAllControlPT = F ## include person-time from controlled trials before end of vaccination refractory period?
    , RCTendOption = 2        ## order to vaccinate unvaccinated invididuals when an RCT ends, see EndTrialFuns.R
    , instVaccDelay = 7 ## delay til instant vacc of everyone after trial ends in trials where delayUnit=0 otherwise
    , small=F ## do a small trial for illustration
    ){
    if(small) {
        numClus <- 4
        clusSize <- 4
    }
    if(maxInfectDay < delayUnit*numClus) stop('maxInfectDay too short. Need enough time to rollout vaccines to all clusters')
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
setHazs <- function(parms=makePop(), browse=F) within(parms, {
    if(browse) browser()
    baseClusHaz <- reParmRgamma(numClus, mean = mu, var = varClus) ## gamma distributed baseline hazards
    dailyDecayRates <- rlnorm(numClus, meanlog = log(weeklyDecay^(1/7)), weeklyDecayVar) 
    daySeq <- seq(-hazIntUnit*ceiling(reordLag/hazIntUnit),maxInfectDay,by=hazIntUnit)
    hazT <- data.table(day = rep(daySeq, each = numClus), cluster = rep(1:numClus, length(daySeq)), clusHaz = 0)
    cHind <- which(names(hazT)=='clusHaz')
    ## mean cluster hazard trajectory
    for(ii in 1:numClus) hazT[which(hazT[,cluster]==ii), clusHaz := baseClusHaz[ii]*dailyDecayRates[ii]^day]
    ## give every individual a lognormally distributed relative risk
    pop$indivRR <- rlnorm(numClus*clusSize, meanlog = 0, sdlog = sdLogIndiv)
    ## create popH which has weekly hazards for all individuals
    popH <- pop[rep(1:nrow(pop), length(daySeq))]
    popH[, day := rep(daySeq, each=nrow(pop))]
    for(dd in daySeq) for(ii in 1:numClus) popH[day==dd & cluster==ii, clusHaz := hazT[day==dd & cluster==ii, clusHaz]]
    popH[, indivHaz := clusHaz*indivRR]
    daySeqLong <- seq(0,maxInfectDay+1000,by=hazIntUnit) ## to avoid problems later
    popHearly <- copy(popH)
    popH <- popH[day >= 0]
    rm(ii, cHind, baseClusHaz, dd, hazT)
})
## setHazs(makePop(makeParms(weeklyDecay=1, weeklyDecayVar=0)))$popH[cluster==1,]
## setHazs(makePop(makeParms(weeklyDecay=.9, weeklyDecayVar=0)))$popH[cluster==1,]

## reparameterize a gamma by mean/var to simulate spatial variation in underlying hazards (change
## later to something more reasonable or based on real data)
reParmRgamma <- function(n, mean, var) {
    if(var > 0) {
        theta <- var/mean
        k <- mean/theta
        rgamma(n, shape = k, scale = theta)
    }else{ ## no variance
        rep(mean, n)
    }
}

reordPop <- function(parms) { ## wrapper around other functions below
    reordFXN <- get(paste0('reord',parms$trial))
    parms <- reordFXN(parms)
    within(parms, { ## return parms
        if(parms$ord!='none') { ## but first, if reordered
            popH[, cluster:=clusIncRank[popH[, cluster]]]
            popH <- arrange(popH, day, cluster)
            popH[,indiv:= rep(1:(numClus*clusSize), length(daySeq))] ## reset indices so they're ordered again by vaccination sequence
            rm(clusIncRank)
        }
        popH$pair <-  NA ## pairs matched for randomization (if matching)
        if(trial=='CRCT' & ord!='none') { ## only paired clusters exist in matched CRCT
            popH$pair <- popH[, cluster %% (numClus/2)]
            popH[pair==0, pair:=numClus/2]
        }
    })
}

reordSWCT <- reordFRCT <- reordRCT <- function(parms) within(parms, {
    if(ord=='BL') { ## if vaccinating highest incidence clusters first (i.e. *NOT* SW randomization of vaccination sequence)
        clusIncRank <- popH[idByClus==1 & day == 0,order(rev(order(clusHaz)))]
    }
    if(ord=='TU') { ## time-updated sequencing, each week vaccinate highest incidence cluster from reordLag days ago
        clusIncRank <- NULL
        for(ii in 1:numClus) { ## for each vaccination day
            dd <- daySeq[ii]
            updatingOrder <- 1:numClus > ii-1 ## i.e. on 3rd day of vaccination, only updating the 3rd vaccination sequence
            currentRank <- popHearly[idByClus==1 & day == dd - reordLag, rev(order(clusHaz))] ## current cluster hazard ordering
            currentRank <- currentRank[!currentRank %in% clusIncRank] ## remove clusters already vaccinated
            clusIncRank <- c(clusIncRank, currentRank[1])
        }
        clusIncRank <- order(clusIncRank)
        rm(currentRank,updatingOrder,dd,ii)
    }
})

reordCRCT <- function(parms) within(parms, {
    if(ord=='BL') { ## if matching clusters on incidence, then randomizing within pairs, then vaccinating highest incidence first
        if(numClus %% 2 == 1) stop("need even # of clusters")
        ## order clusters by mean hazard, in pairs (each row)
        clusIncRank <- matrix(rev(order(popH[day==0 & idByClus==1, clusHaz])), nc = 2, byrow=T) 
        whichVacc <- sample(1:2, nrow(clusIncRank), replace=T) ## which of each pair to vaccinate
        for(ii in 1:nrow(clusIncRank)) if(whichVacc[ii]==2) clusIncRank[ii,1:2] <- clusIncRank[ii,2:1] ## make first column vaccinated group
        clusIncRank <- c(clusIncRank[,1], clusIncRank[,2]) ## reorder of clusters
        clusIncRank <- order(clusIncRank)
        rm(whichVacc, ii)
    }
    ## time-updated sequencing, each week select the two highest incidence pairs that have yet to be
    ## randomized, and randomize one of them to vaccination
    if(ord=='TU') { 
        if(numClus %% 2 == 1) stop("need even # of clusters")
        numPairs <- numClus/2
        clusIncRank <- NULL
        for(ii in 1:min(numPairs,length(daySeq))) { ## for each vaccination day
            dd <- daySeq[ii]
            notRandomized <- (1:numClus)[! 1:numClus %in% as.vector(clusIncRank)] ## haven't already been randomized
            currIncOrd <- notRandomized[rev(order(popHearly[day==dd - reordLag & idByClus==1 & cluster %in% notRandomized, clusHaz]))]
            currentRank <- matrix(currIncOrd, nc = 2, byrow=T)[1,] ## pick top row of paired matrix
            if(rbinom(1,1,.5)) currentRank <- rev(currentRank) ## deterine which of each pair to vaccinate (1st column)
            clusIncRank <- rbind(clusIncRank, currentRank) ## add pair
        }
        clusIncRank <- as.numeric(c(clusIncRank[,1], clusIncRank[,2])) ## vaccinated, control
        clusIncRank <- order(clusIncRank) ## get reording for next line
        rm(notRandomized,dd,ii,numPairs,currIncOrd,currentRank)
    }
})

## p1 <- setHazs(makePop(makeParms('CRCT', clusSize=2, weeklyDecay=.9, weeklyDecayVar=.3, ord='BL')))
## p1 <- reordPop(p1)
## p1$popH
## p1$popH[idByClus==1 & day==0 & cluster <=numClus/2, order(clusHaz)]

setVaccDays <- function(parms) { ## wrapper around other functions below
    setVaccFXN <- get(paste0('set',parms$trial,'vaccDays'))
    parms <- setVaccFXN(parms)
    within(parms, {
        popH$immuneDay <- popH[,vaccDay] + immunoDelay ## vaccine refractory period
        popH$vacc <- popH[, day>=vaccDay]
        popH$immune <- popH[, day>=immuneDay]
        ## reset pop to refrence data table after reordering and then assignment of vaccday stuff
        pop <- select(popH[day==0], indiv, cluster, pair, idByClus, indivRR, vaccDay, immuneDay) 
    })
}
setSWCTvaccDays <- function(parms) within(parms, {
    popH$vaccDay <- delayUnit*(popH[,cluster]-1)
})
setFRCTvaccDays <- setRCTvaccDays <- function(parms) within(parms, { ## assuming same speed rollout as SWCT (unless FRCT)
    popH$vaccDay <- Inf ## unvaccinated
    popH[idByClus > clusSize/2 , vaccDay := delayUnit*(cluster-1)] ## half get vaccinated in each cluster, 1 per interval
})
setCRCTvaccDays <- function(parms) within(parms, {
    popH$vaccDay <- Inf
    popH[cluster <= numClus/2, vaccDay := delayUnit*(cluster-1)] ## first half of clusters (1 per pair) get vaccinated in sequence
})
## To check ordering works
## p1 <- setHazs(makePop(makeParms(clusSize=2, weeklyDecay=.9, weeklyDecayVar=.2, ord='TU', trial='SWCT',small=T)))
## setVaccDays(p1)$popH[,list(cluster,clusHaz, day,vacc,immune)]
## p1 <- setVaccDays(p1)
## p1$popH[idByClus==1,list(cluster,clusHaz, day,vacc,immune)]

## Simulate infections. Takes popH for initial simulation, or EVpopH for end trial vaccination version (requires startInf)
simInfection <- function(parms, whichDo='pop', startInfectingDay = 0, ## startInf can be set to endTrialDay
                         browse = F) 
    within(parms, { 
        if(browse) browser()
        tmp <- get(whichDo)
        tmpH <- get(paste0(whichDo,'H'))
        if(startInfectingDay==0) tmp$infectDay <- tmpH$infectDay <- Inf ## otherwise it's already got some infection data in it
        tmpH[infectDay > startInfectingDay, infectDay := Inf] ## redoing post endDay stuff with additional folks vacc
        tmp[infectDay > startInfectingDay, infectDay := Inf] ## redoing post endDay stuff with additional folks vacc        
        for(dd in daySeq[daySeq>=startInfectingDay]) { ## infection day is beginning of each hazard interval + exponential waiting time
            alreadyInfected <- tmpH[infectDay!=Inf, indiv] ## don't reinfect those already infected
            tmpH[day==dd & !indiv %in% alreadyInfected, 
                 infectDay := dd + rexp(length(indiv), rate = indivHaz*ifelse(immune, 1-vaccEff, 1))] 
            tmpH[day==dd & !indiv %in% alreadyInfected & infectDay > dd + hazIntUnit, 
                 infectDay := Inf] ## reset if it goes into next hazard interval
        }
        ## copy infection days to pop, to use in analysis
        indivInfDays <- tmpH[infectDay!=Inf & infectDay > startInfectingDay, list(indiv,infectDay)]
        indivInfDays <- arrange(indivInfDays, indiv)
        tmp[indiv %in% indivInfDays[,indiv], infectDay:= indivInfDays[,infectDay]]
        assign(whichDo, tmp)
        assign(paste0(whichDo,'H'), tmpH)
        rm(tmp, tmpH, dd,alreadyInfected,indivInfDays)
    })

## p1 <- setHazs(makePop(makeParms(clusSize=300, numClus=20, weeklyDecay=.9, weeklyDecayVar=0, ord='BL')))
## p1 <- reordPop(p1)
## p1 <- simInfection(p1)
## head(p1$pop[infectDay!=Inf, list(cluster, immuneDay, infectDay)],100)

## simulate whole trial and return with all parameters used
simTrial <- function(parms=makeParms(), browse = F) {
    if(browse) browser()
    parms <- makePop(parms) ## make population
    parms <- setHazs(parms) ## set hazards
    parms <- reordPop(parms) ## reorder vaccination sequence by incidence (if applicable)
    parms <- setVaccDays(parms) ## set vaccination days
    parms <- simInfection(parms) ## simulate infection
    return(parms)
}
## simTrial(b=F)
## simTrial(makeParms(small=T),F)

subsArgs <- function(parms, fxn) parms[names(parms) %in% names(formals(fxn))] ## get parameters necessary for a fxn
