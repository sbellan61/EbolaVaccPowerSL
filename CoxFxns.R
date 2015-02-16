library(geepack)

doRelabel <- function(parms, csd, nboot=200, doFXN=doCoxME, verbose = 0) {
    if(verbose==3.45) browser()
    modNm <- paste0('relab', sub('do','',as.character(substitute(doFXN))))
    vee <- doFXN(csd, verbose=verbose)
    effMean <- vee['mean']
    effEsts <- as.numeric(NULL)
    for(bb in 1:nboot) {
        ## randomly reorder the vaccination sequence of clusters, null is their order doesn't affect the
        ## effect size
        if(trial %in% c('SWCT','CRCT'))
            parmsB <- within(parms, {
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
            })
        if(trial %in% c('RCT','FRCT'))
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
            })
        parmsB <- makeSurvDat(parmsB, whichDo='pop', br=F)
        parmsB <- makeGEEDat(parmsB, whichDo = 'popH', verbose=verbose)
        parmsB <- activeFXN(parmsB, whichDo = 'st')
        csdB <- censSurvDat(parmsB, parmsB$maxInfectDay)
        ## summTrial(csdB)[[3]]
        tmpEst <- doFXN(csdB, verbose = 1)['mean']
        effEsts <- c(effEsts, as.numeric(tmpEst))
    }
    pval <- 1 - ecdf(effEsts)(effMean) ## for vaccine being more effective than expected by chance alone
    bootVee <- data.frame(mean=effMean, lci = NA, uci = NA, p = pval, mod=modNm, err=sum(is.na(effEsts)))
    return(bootVee)
}

doBoot <- function(csd, nboot=200, doFXN=doCoxME, verbose = 0) {
    if(verbose==3.4) browser()
    modNm <- paste0('boot', sub('do','',as.character(substitute(doFXN))))
    vee <- doFXN(csd, verbose=verbose)
    effMean <- vee['mean']
    effEsts <- as.numeric(NULL)
    for(bb in 1:nboot) {
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
        tmpEst <- doFXN(csdB, verbose = 1)['mean']
        effEsts <- c(effEsts, as.numeric(tmpEst))
    }
    lci <- quantile(effEsts, .025, na.rm=T)
    uci <- quantile(effEsts, .975, na.rm=T)
    if(is.na(effMean)) effMean <- median(effEsts,na.rm=T)
    bootVee <- data.frame(mean=effMean, lci = lci, uci = uci, p = NA, mod=modNm, err=sum(is.na(lci)))
    return(bootVee)
}

doCoxME <- function(csd, verbose=0) { ## take censored survival object and return vacc effectiveness estimates
    if(verbose==3.1) browser()
    mod <- suppressMessages(try(coxme(Surv(startDay, endDay, infected) ~ immuneGrp + (1|cluster), data = csd), silent=T))
    if(!inherits(mod, 'try-error')) {
        vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        pval <- pnorm(mod$coef/sqrt(vcov(mod)), lower.tail = vaccEffEst[1]>0)*2
        vaccEffEst <- signif(c(vaccEffEst, pval), 3)
        if(vcov(mod)==0) vaccEffEst[2:4] <- NA ## if failing to converge on effect estimate (i.e. 0 variance in beta coefficient)
        if(inherits(mod, 'try-error')) {print('blah'); browser()}
        vaccEffEst <- data.frame(t(vaccEffEst), 'coxME')
        names(vaccEffEst) <- c('mean','lci','uci','p','mod')
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='coxME')
    }
    return(vaccEffEst)
}

doCoxPH <- function(csd, verbose=0) {
    if(verbose==3.2) browser()
    mod <- try(coxph(Surv(startDay, endDay, infected) ~ immuneGrp + cluster(cluster),
                     data=csd), silent=T)
    if(!inherits(mod, 'try-error')) {
        vaccEffEst <- 1-summary(mod)$conf.int['immuneGrp',c(1,4:3)] 
        pval <- summary(mod)$coefficients['immuneGrp','Pr(>|z|)']
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)), signif(pval,3), mod='coxph')
        colnames(vaccEffEst) <- c('mean','lci','uci','p','mod')
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='coxph')
    }
    return(vaccEffEst)
}

doBayesCox <- function(csd, browse = F) {
}

## Cluster level data with one observation per time unit (not for RCTs bc two different covariates
## classes within cluster level, would need individual approach for that)
doGEEclusAR1 <- function(clusDat, verbose=0) { 
    if(verbose==3.6) browser()
    if(trial %in% c('SWCT','CRCT')) {
        mod <- try(geeglm(cases ~ immuneGrp + day, offset = log(atRisk), id = cluster, data = clusDat, family = poisson, corstr = "ar1"), 
                   silent=T)
        if(!inherits(mod, 'try-error')) {
            vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std.err','Pr(>|W|)')])
            vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
            names(vaccEffEst) <- c('mean','lci','uci','p')
            vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GEEClusAR1')
        }else{
            vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEClusAR1')
        }
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEClusAR1')
    }
    return(vaccEffEst)
}

doGLMMclus <- function(clusDat, bayes=T, verbose=0) {
    if(verbose==3.7) browser()
    if(!bayes) mod <- try(glmer(cases ~ immuneGrp + day + (1|cluster) + offset(log(atRisk)), data = clusDat, family = poisson), silent = T)
    if(bayes) mod <- try(bglmer(cases ~ immuneGrp + day + (1|cluster) + offset(log(atRisk)), data = clusDat, family = poisson), silent = T)
    if(!inherits(mod, 'try-error')) {
        vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')])
        vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
        names(vaccEffEst) <- c('mean','lci','uci','p')
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GLMMClus')
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GLMMClus')
    }
    return(vaccEffEst)
}


doGEEsurv <- function(csd, verbose=0) {
    if(verbose==3.5) browser()
    mod <- try(geeglm(infected ~ immuneGrp, data = csd, family = binomial(link='cloglog'), id = cluster, corstr = 'exchangeable'),
               silent=T)
    if(!inherits(mod, 'try-error')) {
        vaccRes <- as.numeric(summary(mod)$coef['immuneGrp', c('Estimate','Std.err','Pr(>|W|)')])
        vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
        names(vaccEffEst) <- c('mean','lci','uci','p')
        vaccEffEst <- data.frame(t(signif(vaccEffEst,3)),  mod='GEEIndivCloglog')
    }else{
        vaccEffEst <- data.frame(mean=NA, lci=NA, uci=NA, p=NA,  mod='GEEIndivCloglog')
    }
    return(vaccEffEst)
}
