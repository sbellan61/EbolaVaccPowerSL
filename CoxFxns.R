doCoxME <- function(csd, verbose=0) { ## take censored survival object and return vacc effectiveness estimates
    if(verbose==3.1) browser()
    mod <- suppressWarnings(coxme(Surv(startDay, endDay, infected) ~ immuneGrp + (1|cluster), data = csd))
    vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
    pval <- pnorm(mod$coef/sqrt(vcov(mod)), lower.tail = vaccEffEst[1]>0)*2
    vaccEffEst <- signif(c(vaccEffEst, pval), 3)
    if(vcov(mod)==0) vaccEffEst[2:4] <- NA ## if failing to converge on effect estimate (i.e. 0 variance in beta coefficient)
    if(inherits(mod, 'try-error')) {print('blah'); browser()}
    vaccEffEst <- data.frame(t(vaccEffEst), 'coxME')
    names(vaccEffEst) <- c('mean','lci','uci','p','mod')
    return(vaccEffEst)
}

doCoxPH <- function(csd, verbose=0) {
    if(verbose==3.2) browser()
    if(trial %in% c('SWCT','CRCT'))
       mod <- try(coxph(Surv(startDay, endDay, infected) ~ immuneGrp + cluster(cluster),
                        data=csd), silent=T)
    if(trial %in% c('RCT','FRCT'))
       mod <- try(coxph(Surv(startDay, endDay, infected) ~ immuneGrp ## + cluster(cluster),
                        data=csd), silent=T)
    vaccEffEst <- 1-summary(mod)$conf.int['immuneGrp',c(1,4:3)] 
    pval <- summary(mod)$coefficients['immuneGrp','Pr(>|z|)']
    vaccEffEst <- data.frame(t(signif(vaccEffEst,3)), signif(pval,3), mod='coxph')
    colnames(vaccEffEst) <- c('mean','lci','uci','p','mod')
    return(vaccEffEst)
}

doGlmer <- function(csd, bayes=T, verbose = 0) {## take censored survival object and return vacc effectiveness estimates using bayesian glme
    if(verbose==3.3) browser()
    if(bayes) 
        mod <- bglmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    if(!bayes) 
        mod <- glmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    
    mod <- glmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), 
                 family=binomial(link='cloglog'),  data = csd,
                 start = list(fixef = c(-4, 0)),
                 control = glmerControl(nAGQ0initStep = FALSE))
    
    vaccRes <- summary(mod)$coefficients['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')] 
    vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
    names(vaccEffEst) <- c('mean','lci','uci','P')
    return(signif(vaccEffEst,3))
}

doBayesCox <- function(csd, browse = F) {

}

doGEE <- function(csd, browse=F) {

}
