## Code adapted from https://github.com/ICI3D/ebola-forecast/ (Carl Pearson)
## Steve Bellan 2015

library(boot); library(ggplot2)
load(file='data/cleanSLData.Rdata')

## weird - what's going on w/ logit transform here?  understand we want to let fitting algo
## explore -/+ inf, but are only interested in decay rates on 0, +inf.  That's not what's actually
## happening w/ logit however
param.xform <- function(x) c(x['rate_0'], decay_rate = as.numeric(exp(x['logdecay_rate'])), nbsize = as.numeric(exp(x['lognbsize'])))
param.xform(c(rate_0=30, logdecay_rate=log(0.01), lognbsize=log(01.2)))

## analysis functions
llgenerator <- function(ratefun, logprobfun, verbose=0) return(
    function(params, cases, times) {
        if(verbose>24) browser()
        ps <- param.xform(params)
        ll <- sum(logprobfun(ps, ratefun(ps, times), cases))
        return(ll)
    }
    )
exp_decay <- function(params, times) params["rate_0"]*exp(-params["decay_rate"]*times)

pois_cases <- function(params, rates, cases, verbose=0) {
    if(verbose>25) browser()
    out <- dpois(cases, lambda=rates, log=TRUE)
    trates <<- rates
    ## print(out)
    return(out)
}
nbinom_cases <- function(params, rates, cases, verbose=0) {
    if(verbose>25) browser()
    out <- dnbinom(cases, size=params["nbsize"], mu=rates, log=TRUE)
    trates <<- rates
    ## print(out)
    return(out)
}
exp_pois_ll <- llgenerator(exp_decay, pois_cases)
exp_nbinom_ll <- llgenerator(exp_decay, nbinom_cases)

mod_exp_decay <- function(params, times) params["rate_0"]*exp(- (params["decay_rate"]*times)^params["shape_parameter"])
mod_exp_pois_ll <- llgenerator(mod_exp_decay, pois_cases)

flat_model <- function(params, times) {
    rep(params["rate_0"],length(times))
}
flat_pois_ll <- llgenerator(flat_model, pois_cases)

expcurve <- function(optimresults, date_zero, params.xform = function(x) x) with(optimresults, {
    ps <- params.xform(par)
                                        #print(class(date_zero))
    return(function(d) {
        t <- d - date_zero
        ps["rate_0"]*exp(-ps["decay_rate"]*t)
    })
})

params <- function(src, censor_interval, include_interval, 
                   param.guess = c(rate_0=30, logdecay_rate=log(0.01)), verbose = 0,
                   lognbsize=log(01.2), ll = 'exp_nbinom_ll') {
    if(verbose>22) browser()
    end_date <- max(src$Date) - censor_interval
    start_date <- end_date - include_interval
    slice <- src[(start_date <= src$Date) & (src$Date <= end_date),] ## data.table doesn't work with shiny for some reason?
    slice$t <- as.numeric(slice$Date - min(slice$Date))
    if(sum(slice$int>1)>0) print('not always daily incidence')
    if(ll=='exp_pois_ll') parscale=c(1,1/10)
    if(ll=='exp_nbinom_ll') {
        parscale=c(1,1/10,1/10)
        param.guess <- c(rate_0=30, logdecay_rate=log(0.01), lognbsize = lognbsize)
    }
    llfxn <- get(ll)
    res <- optim(param.guess, llfxn, gr=NULL, cases = slice$cases, times = slice$t, control=list(fnscale=-1, parscale=parscale), hessian=TRUE)
    vcmat <- solve(-res$hessian)
    res$par <- param.xform(res$par)
    cis <- with(res,{
        offset <- 1.96*sqrt(vcmat[2,2])
        par["decay_rate"] + c(low = offset, est=0, hi = -offset) # hi b = hi decay = low time
    })
    list(t0 = start_date, par=res$par, cis=cis)
}

doProj <- function(src, forecast_date = as.Date('2015-12-31'), censor_interval = 0, include_interval = 30, minCases = 15, 
                   beforeDays = 30, moreDays = 90, minDecayRate = .01, ll = 'exp_nbinom_ll', model = expcurve,
                   verbose=0) {
    if(verbose>20) browser()
    endDate <- max(src$Date) - censor_interval
    ## minimum cases in last interval to fit exponential, otherwise increase to max # cases
    fromMax <- T
    if(fromMax) startDate <- src[length(cases) + 1 - which.max(rev(cases)), Date]
    include_interval <- endDate-startDate
    fit_ref <- params(src, censor_interval, include_interval, ll = ll, verbose=verbose)
    fit_ref$par['decay_rate'] <- max(fit_ref$par['decay_rate'], minDecayRate)
    srcProj <- data.table(Date = seq(min(src$Date), forecast_date, by='day'))
    maxDataDate <- max(src$Date)
    src <- merge(src, srcProj, by = 'Date', all=T)
    src[Date < maxDataDate & is.na(cases), cases:=0] ## fill in dates with missing values for plotting
    src[, reg := reg[1]]
    src$days <- as.numeric(src$Date)
    src$fit <- model(fit_ref, date_zero = src[Date==startDate, days])(src$days)
    return(list(fit = fit_ref, src = src[,list(Date,reg, cases,fit)], startDate = startDate, endDate = endDate,
                moreDays=moreDays, beforeDays=beforeDays, include_interval=include_interval, model = model))
}

forecast <- function(fit, main=NULL, nbsize = NULL, doPlot = T, xticks = T,  ylim = NULL, xlim = NULL, verbose=0) with(fit, {
    if(verbose>20) browser()
    if(!is.null(xlim)) {
        startX <- xlim[1]
        endX <- xlim[2] 
    }else{
        startX <- startDate - beforeDays
        endX <- endDate + moreDays    
        xlim <- c(startX, endX)
    }
    if(is.null(nbsize)) nbsize <- fit$par$nbsize ## if not specified in arguments
    src$proj <- src[, rnbinom(length(fit), mu = fit, size  = nbsize)]
    if(doPlot) {
        src[Date < endDate & Date > startX, plot(Date, cases, type = 'h', bty = 'n', las = 2, xaxt = 'n', 
                                                 xlim=xlim, ylim=ylim, xlab = '', lwd = 3)]
        rgDates <- seq.Date(startX, endX, by = 1)
        seqDates <- rgDates[format(rgDates, '%d') %in% c('01')]
        seqDatesMid <- rgDates[format(rgDates, '%d') %in% c('15')]
        axis.Date(1, at = seqDatesMid, labels = NA)
        if(xticks) axis.Date(1, at = seqDates, format = '%b', las = 2, tck=0)
        rect(startDate, 0, endDate, par('usr')[4], col = rgb(0,.5,0,.3), border=NA)
        title(main)
        src[Date > endDate, lines(Date, proj, type = 'h', lwd = 3, col = 'red')]
        src[Date > startDate, lines(Date, fit, lty = 1, col = 'dodger blue', lwd = 2)]
        rect(endDate, 0, endDate + moreDays, par('usr')[4], col = rgb(0.5,0,0,.3), border=NA)
        text(startDate+include_interval/2, par('usr')[4], 'fitting', pos=1)
        text(endDate+moreDays/2, par('usr')[4], 'forecasting', pos=1)
    }
    return(src)
})

createHazTraj <- function(fits, nbsize = 1.2, trialStartDate = as.Date('2015-02-01'), xlim = as.Date(c('2014-09-15','2015-12-01')),
                          propInTrial = .03, numClus = 20, clusSize = 300, weeks = T, verbose=0) {
    hazTList <- NULL
    if(verbose>20) browser()
    for(cc in 1:numClus) {
        fit <- fits[[sample(regs, 1)]]
        src <- forecast(fit, doPlot = F, nbsize = nbsize, xlim = xlim, verbose=verbose)
        src$day <- src[, as.numeric(Date - trialStartDate)]
        lastDataDay <- src[max(which(!is.na(src$cases))), day]
        src[day < lastDataDay & is.na(cases), cases := 0] ## fill in over interval without reporting
        src$haz <- src$cases
        src[day > lastDataDay, haz := proj]
        src[, haz := haz * propInTrial / clusSize]
        src[day > -60, list(day, haz)]
        if(weeks) {
            src$week <- src[, floor(day/7)]
            src <- src[, list(haz = mean(haz, na.rm=T), Date = min(Date)), week] ## get mean daily hazard by week
        }
        src$cluster <- cc
        hazTList[[cc]] <- as.data.frame(src)
    }
    hazT <- rbindlist(hazTList)
    setnames(hazT, 'haz', 'clusHaz')
    hazT$day <- hazT[, week*7]
    hazT <- arrange(hazT, day, cluster)
    return(hazT[, list(cluster, day, Date, clusHaz)])
}

makeTransparent<-function(someColor, alpha=150) {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
