doCoxPH <- function(csd, pkg='coxme', browse=F) { ## take censored survival object and return vacc effectiveness estimates
    if(browse) browser()
    if(pkg=='coxph') {
        mod <- try(coxph(Surv(startDay, endDay, infected) ~ immuneGrp + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
                         outer.max=1000, iter.max=10000,
                         data=csd), silent=T)
        vaccEffEst <- 1-summary(mod)$conf.int['immuneGrp',c(1,4:3)] 
        pval <- summary(mod)$coefficients['immuneGrp','p']
        vaccEffEst <- c(vaccEffEst, pval)
    }
    if(pkg=='coxme') {
        mod <- suppressWarnings(coxme(Surv(startDay, endDay, infected) ~ immuneGrp + (1|cluster), data = csd))
        vaccEffEst <- 1-exp(mod$coef + c(0, 1.96, -1.96)*sqrt(vcov(mod)))
        pval <- pnorm(mod$coef/sqrt(vcov(mod)))*2
        vaccEffEst <- c(vaccEffEst, pval)
        if(vcov(mod)==0) vaccEffEst[2:4] <- NA ## if failing to converge on effect estimate (i.e. 0 variance in beta coefficient)
    }
    ##     mod <- coxph(Surv(startDay, endDay, infected) ~ immuneGrp, data=csd) ## without frailty
    if(inherits(mod, 'try-error')) {print('blah'); browser()}
    names(vaccEffEst) <- c('mean','lci','uci','p')
    return(signif(vaccEffEst,3))
}

doGlmer <- function(csd, bayes=F, browse = F) {## take censored survival object and return vacc effectiveness estimates using bayesian glme
    if(browse) browser()
    if(bayes) mod <- bglmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    if(!bayes) mod <- glmer(infected ~ immuneGrp + (1|cluster) + offset(log(perstime)), family=binomial(link='cloglog'),  data = csd)
    vaccRes <- summary(mod)$coefficients['immuneGrp', c('Estimate','Std. Error','Pr(>|z|)')] 
    vaccEffEst <- c(1 - exp(vaccRes[1] + c(0, 1.96, -1.96) * vaccRes[2]), vaccRes[3])
    names(vaccEffEst) <- c('mean','lci','uci','P')
    return(signif(vaccEffEst,3))
}

doBayesCox <- function(csd, browse = F) {

}
