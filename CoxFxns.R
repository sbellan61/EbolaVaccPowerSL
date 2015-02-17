library(geepack)

## make the highest hazard uninfected individuals in a vacc and control group be infected to allow
## for conservative CI calculation when 0 infections in certain groups
infBump <- function(parms, verbose = 0) {
    if(verbose==3.9) browser()
    if(verbose>0) print('bumping infections to deal with divergence')
    parmsE <- copy(parms)
    parmsE <- within(parmsE, {
        infecteds <- popH[, list(infected = sum(infectDay!=Inf)), indiv]
        uninfecteds <- infecteds[infected==0, indiv]
        newInfs <- popH[indiv %in% uninfecteds,
                        list(indiv = indiv[indivHaz==max(indivHaz)], 
                             day = day[indivHaz==max(indivHaz)]), immune]
        popH[indiv==newInfs[1,indiv] & day==newInfs[1,day], infectDay := day + 6.9]
        popH[indiv==newInfs[2,indiv] & day==newInfs[2,day], infectDay := day + 6.9]
        indivInfDays <- popH[infectDay!=Inf & indiv %in% newInfs[,indiv], list(indiv,infectDay)]
        indivInfDays <- arrange(indivInfDays, indiv)
        pop[indiv %in% indivInfDays[,indiv], infectDay:= indivInfDays[,infectDay]]
        pop[indiv %in% indivInfDays[,indiv], ]
        rm(infecteds, uninfecteds, newInfs, indivInfDays)
    })
    parmsE <- makeSurvDat(parmsE, whichDo='pop', br=F)
    parmsE <- makeGEEDat(parmsE, whichDo = 'popH', verbose=verbose)
    parmsE <- activeFXN(parmsE, whichDo = 'st')
    return(parmsE)
}

doRelabel <- function(parms, csd, bump=F, nboot=200, doFXN=doCoxME, verbose = 0, verbFreqRelab=10) {
    if(verbose==3.45) browser()
    if(verbose>0) print('doing relabeled models')
    modNm <- paste0('relab', sub('do','',as.character(substitute(doFXN))))
    vee <- doFXN(csd, verbose=0)
    effMean <- vee['mean']
    effEsts <- as.numeric(NULL)
    for(bb in 1:nboot) {
        if(verbose>.5 & (bb %% verbFreqRelab == 0)) print(paste('on',bb,'of',nboot))
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
                popH[, vacc := day>= vaccDay]
                popH[, immune := day>= immuneDay]
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
        parmsB <- makeSurvDat(parmsB, whichDo='pop', br=F)
        parmsB <- makeGEEDat(parmsB, whichDo = 'popH', verbose=verbose)
        parmsB <- activeFXN(parmsB, whichDo = 'st')
        csdB <- censSurvDat(parmsB, parmsB$maxInfectDay)
        ## summTrial(csdB)[[3]]
        tmpEst <- doFXN(csdB, verbose = 0)['mean']
        effEsts <- c(effEsts, as.numeric(tmpEst))
    }
    pval <- ecdf(effEsts)(effMean) ## for vaccine being more effective than expected by chance alone
    if(effMean > 0) pval <- 1 - pval
    bootVee <- data.frame(mean=effMean, lci = NA, uci = NA, p = pval, mod=modNm, bump = NA, err=sum(is.na(effEsts)))
    bootVee <- bumpAdjust(bootVee, csd, bump, nonpar=T)
    return(bootVee)
}

doBoot <- function(csd, nboot=200, bump=F, doFXN=doCoxME, verbose = 0, verbFreqBoot=10) {
    if(verbose==3.4) browser()
    if(verbose>0) print('bootstrapping')
    modNm <- paste0('boot', sub('do','',as.character(substitute(doFXN))))
    vee <- doFXN(csd, verbose=0)
    effMean <- vee['mean']
    effEsts <- as.numeric(NULL)
    for(bb in 1:nboot) {
        if(verbose>.5 & (bb %% verbFreqBoot == 0)) print(paste('on',bb,'of',nboot))
        if(trial=='CRCT' & ord!='none')
            bootby <- csd[,unique(pair)]    else    bootby <- csd[,unique(cluster)]
        clsB <- sample(bootby, length(bootby), replace=T)
        clsB <- clsB[order(clsB)]
        clsB <- table(factor(clsB, bootby))
        csdB <- copy(csd)
        csdB$reps <- NA
        if(trial=='CRCT' & ord!='none')
            csdB[, reps := clsB[pair]] else csdB[, reps := clsB[cluster]] 
        csdB <- csdB[rep(1:length(reps), reps)]
        tmpEst <- doFXN(csdB, verbose = 0)['mean']
        effEsts <- c(effEsts, as.numeric(tmpEst))
    }
    lci <- quantile(effEsts, .025, na.rm=T)
    uci <- quantile(effEsts, .975, na.rm=T)
    if(is.na(effMean)) effMean <- median(effEsts,na.rm=T)
    bootVee <- data.frame(mean=effMean, lci = lci, uci = uci, p = NA, mod=modNm, bump= NA, err=sum(is.na(lci)))
    bootVee <- bumpAdjust(bootVee, csd, bump, nonpar=T)
    return(bootVee)
}

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
            if(!nonpar) vee['p'] <- NA
        }
    }
    return(vee)
}

doCoxME <- function(csd, bump = F, verbose=0) { ## take censored survival object and return vacc effectiveness estimates
    if(verbose==3.1) browser()
    if(verbose>0) print('fitting vanilla coxME')
    mod <- suppressMessages(try(coxme(Surv(startDay, endDay, infected) ~ immuneGrp + (1|cluster), data = csd), silent=T))
    if(!inherits(mod, 'try-error')) {
        vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        pval <- pnorm(mod$coef/sqrt(vcov(mod)), lower.tail = vaccEffEst[1]>0)*2
        vaccEffEst <- signif(c(vaccEffEst, pval), 3)
        if(vcov(mod)==0) vaccEffEst[2:4] <- NA ## if failing to converge on effect estimate (i.e. 0 variance in beta coefficient)
        vaccEffEst <- data.frame(t(vaccEffEst), 'coxME', bump)
        names(vaccEffEst) <- c('mean','lci','uci','p','mod', 'bump')
        vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='coxME', bump = bump)
    }
    return(vaccEffEst)
}

## Cluster level data with one observation per time unit (not for RCTs bc two different covariates
## classes within cluster level, would need individual approach for that)
doGEEclusAR1 <- function(clusDat, csd, bump=F, verbose=0) { 
    if(verbose==3.6) browser()
    if(verbose>0) print('fitting GEEclusAR1')
    if(trial %in% c('SWCT','CRCT')) {
        mod <- try(geeglm(cases ~ immuneGrp + day, offset = log(atRisk), id = cluster, data = clusDat, family = poisson, corstr = "ar1"), 
                   silent=T)
        if(!inherits(mod, 'try-error')) {
            vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std.err','Pr(>|W|)')])
            vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
            names(vaccEffEst) <- c('mean','lci','uci','p')
            vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GEEClusAR1', bump=bump)
            vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
        }else{
            vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEClusAR1')
        }
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEClusAR1')
    }
    return(vaccEffEst)
}

doGLMMclus <- function(clusDat, csd, bump=F, bayes=T, verbose=0) {
    if(verbose==3.7) browser()
    if(verbose>0) print('fitting GLMMclus')
    if(!bayes) mod <- try(glmer(cases ~ immuneGrp + day + (1|cluster) + offset(log(atRisk)), data = clusDat, family = poisson), silent = T)
    if(bayes) mod <- try(bglmer(cases ~ immuneGrp + day + (1|cluster) + offset(log(atRisk)), data = clusDat, family = poisson), silent = T)
    if(!inherits(mod, 'try-error')) {
        vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')])
        vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
        names(vaccEffEst) <- c('mean','lci','uci','p')
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GLMMClus', bump = bump)
            vaccEffEst <- bumpAdjust(vaccEffEst, csd, bump)
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GLMMClus')
    }
    return(vaccEffEst)
}
