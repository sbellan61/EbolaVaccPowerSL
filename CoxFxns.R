####################################################################################################
## Functions to analyze simulated trial data.
####################################################################################################
## Code base accompanying:
## 
## Bellan, SE, JRC Pulliam, CAB Pearson, DChampredon, SJ Fox, L Skrip, AP Galvani, M Gambhir, BA
## Lopman, TC Porco, LA Meyers, J Dushoff (2015). The statistical power and validity of Ebola
## vaccine trials in Sierra Leone: A simulation study of trial design and analysis. _Lancet
## Infectious Diseases_.
##
## Steve Bellan, March 2015
## License at bottom.
####################################################################################################

library(geepack)

## make the highest hazard uninfected individuals in a vacc and control group be infected to allow
## for conservative CI calculation when 0 infections in certain groups
infBump <- function(parms) {
    if(parms$verbose==3.9) browser()
    if(parms$verbose>0) print('bumping infections to deal with divergence')
    parmsE <- copy(parms)
    parmsE <- within(parmsE, {
        infecteds <- popH[, list(infected = sum(infectDay!=Inf)), indiv]
        uninfecteds <- infecteds[infected==0, indiv]
        popH$immuneDayThink <- popH[, vaccDay + immunoDelay]
        popH$firstActive <- 0
        if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
            popH[, firstActive := min(immuneDayThink), by = pair]
        if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
            popH[, firstActive := min(immuneDayThink), cluster]
        popH$active <- popH[,day>=firstActive]
        newInfs <- arrange(popH[indiv %in% uninfecteds & active==T], desc(indivHaz))
        contInf <- newInfs[immune==F, list(indiv = indiv[indivHaz==max(indivHaz)], day = day[indivHaz==max(indivHaz)])][1,]
        vaccInf <- newInfs[immune==T & indiv!=contInf[,indiv], 
                           list(indiv = indiv[indivHaz==max(indivHaz)], day = day[indivHaz==max(indivHaz)])][1,]
        newInfs <- rbind(contInf,vaccInf)
        popH[indiv==newInfs[1,indiv] & day==newInfs[1,day], infectDay := day + 6.9]
        popH[indiv==newInfs[2,indiv] & day==newInfs[2,day], infectDay := day + 6.9]
        indivInfDays <- popH[infectDay!=Inf & indiv %in% newInfs[,indiv], list(indiv,infectDay)]
        indivInfDays <- arrange(indivInfDays, indiv)
        pop[indiv %in% indivInfDays[,indiv], infectDay:= indivInfDays[,infectDay]]
        pop[indiv %in% indivInfDays[,indiv], ]
        rm(infecteds, uninfecteds, newInfs, indivInfDays, st, stActive, clusDat)
    })
    parmsE <- makeSurvDat(parmsE)
    parmsE <- makeGEEDat(parmsE)
    parmsE <- activeFXN(parmsE)
    return(parmsE)
}

##  permutation P value
permP <- function(x,na.rm=T) min(mean(x>=x[1],na.rm=na.rm), mean(x<=x[1],na.rm=na.rm))

##  permutation tests
doRelabel <- function(parms, csd, bump=F, nboot=200, doMods=modsToDo, verbFreqRelab=10, minCases=0) {
    if(parms$verbose==3.45) browser()
    if(parms$verbose>0) print(paste('doing relabeled models:', paste(unlist(modsToDo), collapse=', ')))
    doFXNs <- lapply(doMods, function(x) get(paste0('do',x)))
    names(doFXNs) <- doMods
    nmods <- length(doMods)
    veePerm <- data.table(matrix(0,nboot,nmods))
    setnames(veePerm, unlist(doMods))
    ## First row in table is effect estimate from simulation, remainin nboot-1 are from permutations
    for(ii in 1:nmods) set(veePerm, i=1L, j=ii, doFXNs[[ii]](parms=within(parms, {verbose=0}), csd=csd, bump=F)[1,'mean'])
    numCases <- csd[,sum(infectDay!=Inf)]
    if(numCases >= minCases) { #
        for(bb in 2:nboot) {
            if(parms$verbose>.5 & (bb %% verbFreqRelab == 0)) print(paste('on',bb,'of',nboot))
            ## randomly reorder the vaccination sequence of clusters, null is their order doesn't affect the
            ## effect size
            if(trial %in% c('SWCT','CRCT')) {
                parmsB <- copy(parms)
                parmsB <- within(parmsB, {
                    relabDT <- pop[!duplicated(cluster), list(cluster, vaccDay)]
                    relabDT[, vaccDay := sample(vaccDay, size = length(vaccDay), replace = F)]
                    pnms <- colnames(pop)
                    pop <- merge(pop[,!'vaccDay',with=F], relabDT, by = 'cluster')
                    pop[, immuneDay:=vaccDay+immunoDelay]
                    pop <- setcolorder(pop, pnms)
                    pnms <- colnames(popH)
                    popH <- merge(popH[,!'vaccDay',with=F], relabDT, by = 'cluster')
                    popH[, immuneDay := vaccDay + immunoDelay]
                    popH[, vacc := day >= vaccDay]
                    popH[, immune := day >= immuneDay]
                    popH <- setcolorder(popH, pnms)
                    rm(relabDT, pnms) 
                }) }
            if(trial %in% c('RCT','FRCT')) {
                parmsB <- copy(parms)
                parmsB <- within(parms, {
                    ## Relabel who got vaccinated within a cluster
                    whoVacc <- sample(1:clusSize, clusSize/2, replace=F)
                    pop[idByClus %in% whoVacc, vaccDay := min(vaccDay), cluster] 
                    pop[!idByClus %in% whoVacc, vaccDay := Inf]
                    pop[, immuneDay:=vaccDay+immunoDelay]
                    pnms <- colnames(popH)
                    popH <- merge(popH[,!'vaccDay',with=F], pop[,list(indiv,vaccDay)], by = 'indiv')
                    popH <- arrange(popH, day, indiv)
                    popH[, immuneDay := vaccDay + immunoDelay]
                    popH[, vacc := day>= vaccDay]
                    popH[, immune := day>= immuneDay]
                    popH <- setcolorder(popH, pnms)
                    rm(relabDT, pnms)
                }) }
            parmsB <- makeSurvDat(parmsB, whichDo='pop')
            parmsB <- makeGEEDat(parmsB, whichDo = 'popH')
            parmsB <- activeFXN(parmsB, whichDo = 'st')
            csdB <- censSurvDat(parmsB)
            parmsB$verbose <- 0 ## don't want printouts within resampling
            for(ii in 1:nmods) set(veePerm, i=as.integer(bb), j=ii, doFXNs[[ii]](parms=parmsB, csd=csdB, bump=F)[1,'mean'])
        }
        bootVee <- data.frame(mean = as.numeric(veePerm[1]), lci = NA, uci = NA, p = apply(veePerm, 2, permP), 
                              mod = paste0('relab',unlist(doMods)), bump = F, err = colSums(is.na(veePerm)))
    }else{ ## not enough cases
        bootVee <- data.frame(mean = as.numeric(veePerm[1]), lci = NA, uci = NA, p = NA,
                              mod = paste0('relab',unlist(doMods)), bump = F, err = NA)
    }
    return(bootVee)
}

##  cluster-level bootstrap
doBoot <- function(parms, csd, nboot=200, bump=F, doMods=modsToDo, verbFreqBoot=10, minCases=0) {
    if(parms$verbose==3.4) browser()
    if(parms$verbose>0) print(paste('doing bootstrapped models:', paste(unlist(modsToDo), collapse=', ')))
    doFXNs <- lapply(doMods, function(x) get(paste0('do',x)))
    names(doFXNs) <- doMods
    nmods <- length(doMods)
    veeBoot <- data.table(matrix(0,nboot,nmods))
    setnames(veeBoot, unlist(doMods))
    ## First row in table is effect estimate from simulation, remainin nboot-1 are from permutations
    for(ii in 1:nmods) set(veeBoot, i=1L, j=ii, doFXNs[[ii]](parms=within(parms, {verbose=0}), csd=csd, bump=F)[1,'mean'])
    numCases <- csd[,sum(infectDay!=Inf)]
    if(numCases >= minCases) {
        for(bb in 2:nboot) {
            if(parms$verbose>.5 & (bb %% verbFreqBoot == 0)) print(paste('on',bb,'of',nboot))
            if(trial=='CRCT' & ord!='none')
                bootby <- csd[,unique(pair)]    else    bootby <- csd[,unique(cluster)]
            clsB <- sample(bootby, length(bootby), replace=T)
            clsB <- clsB[order(clsB)]
            clsB <- table(factor(clsB, bootby))
            csdB <- copy(csd)
            parmsB <- copy(parms)
            parmsB$clusDat$reps <- csdB$reps <- NA
            if(trial=='CRCT' & ord!='none') {
                csdB[, reps := clsB[pair]] 
                parmsB$clusDat[, reps := clsB[pair]] 
            }else{ 
                csdB[, reps := clsB[cluster]] 
                parmsB$clusDat[, reps := clsB[cluster]] 
            }
            csdB <- csdB[rep(1:length(reps), reps)]
            parmsB$clusDat <- parmsB$clusDat[rep(1:length(reps), reps)]
            parmsB$verbose <- 0 ## don't want printouts within resampling
            for(ii in 1:nmods) set(veeBoot, i=as.integer(bb), j=ii, doFXNs[[ii]](parms=parmsB, csd=csdB, bump=F)[1,'mean'])
        }
        bootVee <- data.frame(mean = as.numeric(veeBoot[1])
                              , lci = apply(veeBoot[-1], 2, function(x) quantile(x,.025,na.rm=T)) ## [-1] exclude real estimate
                              , uci = apply(veeBoot[-1], 2, function(x) quantile(x,.975,na.rm=T)) 
                              , p = NA
                              , mod = paste0('boot',unlist(doMods)), bump = F, err = colSums(is.na(veeBoot)))
    }else{
        bootVee <- data.frame(mean = as.numeric(veeBoot[1])
                              , lci = NA
                              , uci = NA
                              , p = NA
                              , mod = paste0('boot',unlist(doMods)), bump = F, err = NA)
    }
    return(bootVee)
}

##  If running infection bumps,  adjust confidence interval bounds accordingly.
bumpAdjust <- function(vee, csd, bump, nonpar=F) {
    if(bump) {
        zeroVacc <- csd[,list(oneBumped = sum(infectDay!=Inf)==1), immuneGrp]
        if(zeroVacc[immuneGrp==1, oneBumped]) { ## there were zero cases in the vaccinated group
            vee['mean'] <- 1
            if(!nonpar) vee['uci'] <- 1
        }
        if(zeroVacc[immuneGrp==0, oneBumped]) { ## there were zero cases in the vaccinated group
            vee['mean'] <- -Inf
            if(!nonpar) vee['lci'] <- -Inf
        }
        if(sum(zeroVacc[,oneBumped==1])==2) { ## there were zero cases in both groups
            vee['mean'] <- NA
            if(!nonpar) vee['lci'] <- NA
            if(!nonpar) vee['uci'] <- NA
            if(!nonpar) vee['p'] <- 1
        }
    }
    return(vee)
}

## Cox PH gamma frailty
doCoxME <- function(parms, csd, bump = F) { ## take censored survival object and return vacc effectiveness estimates
    if(parms$verbose==3.1) browser()
    if(parms$verbose>0) print('fitting vanilla coxME')
    mod <- try(coxme(Surv(startDay, endDay, infected) ~ immuneGrp + (1|cluster), data = csd), silent=T)
    if(!inherits(mod, 'try-error')) {
        vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        pval <- pnorm(mod$coef/sqrt(vcov(mod)), lower.tail = vaccEffEst[1]>0)*2
        vaccEffEst <- signif(c(vaccEffEst, pval), 3)
        if(is.na(vcov(mod)) | vcov(mod)==0) 
            vaccEffEst[2:4] <- NA ## if failing to converge on effect estimate (i.e. 0 variance in beta coefficient)
        vaccEffEst <- data.frame(t(vaccEffEst), 'coxME', bump = bump, err = 0)
        names(vaccEffEst) <- c('mean','lci','uci','p','mod', 'bump', 'err')
        vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='coxME', bump = bump, err = 1)
    }
    return(vaccEffEst)
}

## Cluster level data with one observation per time unit (not for RCTs bc two different covariates
## classes within cluster level, would need individual approach for that)
doGEEclusAR1 <- function(parms, csd, bump=F) { ## PROBLEMATIC, crashed on cluster repeatedly, failed to converge, exhibited poor coverage & inflated false positives
    if(parms$verbose==3.6) browser()
    if(parms$verbose>0) print('fitting GEEclusAR1')
    if(trial %in% c('SWCT','CRCT')) {
        mod <- try(geeglm(cases ~ immuneGrp + day, offset = log(atRisk), id = cluster, data = parms$clusDat, 
                          family = poisson, corstr = "ar1"), 
                   silent=T)
        if(!inherits(mod, 'try-error')) {
            vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std.err','Pr(>|W|)')])
            vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
            names(vaccEffEst) <- c('mean','lci','uci','p')
            vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GEEClusAR1', bump=bump, err = 0)
            vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
        }else{
            vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEClusAR1', bump = bump, err = 1)
        }
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEClusAR1', bump= bump, err = NA)
    }
    return(vaccEffEst)
}

##  Mixed effects Poisson regression
doGLMMclusFr <- function(parms, csd, bump=F) doGLMMclus(parms, csd, bump, bayes=F) 
doGLMMclusBy <- function(parms, csd, bump=F) doGLMMclus(parms, csd, bump, bayes=T) 
doGLMMclus <- function(parms, csd, bump=F, bayes=T) {
    if(parms$verbose==3.7) browser() ## PROBLEMATIC, crashed on cluster repeatedly (segmentation faults), failed to converge
    if(parms$verbose>0) print('fitting GLMMclus')
    if(!bayes) mod <- try(glmer(cases ~ immuneGrp + day + (1|cluster) + offset(log(atRisk)), 
                                data = parms$clusDat, family = poisson), silent = T)
    if(bayes) mod <- try(bglmer(cases ~ immuneGrp + day + (1|cluster) + offset(log(atRisk)), 
                                data = parms$clusDat, family = poisson), silent = T)
    if(!inherits(mod, 'try-error')) {
        vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')])
        vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
        names(vaccEffEst) <- c('mean','lci','uci','p')
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GLMMclus', bump = bump, err = 0)
        vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GLMMclus', bump = NA, err = 1)
    }
    return(vaccEffEst)
}

##  Fixed effects Poisson regression
doGLMclus <- function(parms, csd, bump=F) {
    if(parms$verbose==3.75) browser()
    if(parms$verbose>0) print('fitting GLMclus')
    mod <- try(glm(cases ~ immuneGrp + day + offset(log(atRisk)), data = parms$clusDat, family = poisson), silent = T)
    if(!inherits(mod, 'try-error')) {
        vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')])
        vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
        names(vaccEffEst) <- c('mean','lci','uci','p')
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GLMclus', bump = bump, err = 0)
        vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GLMclus', bump = NA, err = 1)
    }
    return(vaccEffEst)
}

##  Fixed effects Poisson regression with cluster-level fixed effects
doGLMFclus <- function(parms, csd, bump=F) {
    if(parms$verbose==3.76) browser()
    if(parms$verbose>0) print('fitting GLMFclus')
    mod <- try(glm(cases ~ immuneGrp + day + factor(cluster) + offset(log(atRisk)), 
                   data = parms$clusDat, family = poisson), silent = T)
    if(!inherits(mod, 'try-error')) {
        vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')])
        vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
        names(vaccEffEst) <- c('mean','lci','uci','p')
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GLMFclus', bump = bump, err = 0)
        vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GLMFclus', bump = NA, err = 1)
    }
    return(vaccEffEst)
}

####################################################################################################
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/sbellan61/EbolaVaccPowerSL
####################################################################################################
