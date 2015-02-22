## Code adapted from https://github.com/ICI3D/ebola-forecast/ (Carl Pearson)
## Steve Bellan 2015

library(boot); library(ggplot2)
load(file='data/cleanSLData.Rdata')

## weird - what's going on w/ logit transform here?  understand we want to let fitting algo
## explore -/+ inf, but are only interested in decay rates on 0, +inf.  That's not what's actually
## happening w/ logit however
param.xform <- function(x) c(x['rate_0'], decay_rate = as.numeric(exp(x['logdecay_rate'])), nbsize = as.numeric(exp(x['lognbsize'])))
param.xform(c(rate_0=30, logdecay_rate=log(0.01), lognbsize=log(0.9)))

## analysis functions
llgenerator <- function(ratefun, logprobfun) return(
    function(params, cases, times) {
        ps <- param.xform(params)
        ll <- sum(logprobfun(ps, ratefun(ps, times), cases))
        return(ll)
    }
)
exp_decay <- function(params, times) params["rate_0"]*exp(-params["decay_rate"]*times)

pois_cases <- function(params, rates, cases) {
    out <- dpois(cases, lambda=rates, log=TRUE)
    trates <<- rates
    ## print(out)
    return(out)
}
nbinom_cases <- function(params, rates, cases) {
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

params <- function(src, censor_interval, include_interval, param.guess = c(rate_0=30, logdecay_rate=log(0.01), lognbsize=log(0.9)), browse = F, ll = exp_nbinom_ll) {
    if(browse) browser()
    end_date <- max(src$Date) - censor_interval
    start_date <- end_date - include_interval
    slice <- src[(start_date <= src$Date) & (src$Date <= end_date),] ## data.table doesn't work with shiny for some reason?
    slice$t <- as.numeric(slice$Date - min(slice$Date))
    if(sum(slice$int>1)>0) print('not always daily incidence')
    res <- optim(param.guess, ll, gr=NULL, cases = slice$cases, times = slice$t, control=list(fnscale=-1, parscale=c(1,1/10,1/10)), hessian=TRUE)
    vcmat <- solve(-res$hessian)
    res$par <- param.xform(res$par)
    cis <- with(res,{
        offset <- 1.96*sqrt(vcmat[2,2])
        par["decay_rate"] + c(low = offset, est=0, hi = -offset) # hi b = hi decay = low time
    })
    list(t0 = start_date, par=res$par, cis=cis)
}

doProj <- function(src, forecast_date = tail(src$Date,1), censor_interval = 0, include_interval = 30, minCases = 15, 
                      beforeDays = 30, moreDays = 90, minDecayRate = .01, ll = exp_nbinom_ll, model = expcurve) {
    endDate <- max(src$Date) - censor_interval
    # startDate <- endDate - include_interval
    # startX <- startDate - beforeDays
    ## minimum cases in last interval to fit exponential, otherwise increase to max # cases
    fromMax <- T #ifelse(src[Date > max(Date)-include_interval, sum(cases)] < minCases, T, F)
    if(fromMax) startDate <- src[length(cases) + 1 - which.max(rev(cases)), Date]
    include_interval <- endDate-startDate
    fit_ref <- params(src, censor_interval, include_interval, ll = ll)
    fit_ref$par['decay_rate'] <- max(fit_ref$par['decay_rate'], minDecayRate)
    src$fit <- model(fit_ref, date_zero = as.numeric(startDate))(as.numeric(src$Date))
    return(list(fit = fit_ref, src = src[,list(Date,reg, cases,fit)], startDate = startDate, endDate = endDate,
                moreDays=moreDays, beforeDays=beforeDays, include_interval=include_interval, model = model))
}

forecast <- function(fit, main=NULL, nbsize = .9, doPlot = T, xlim = NULL) with(fit, {
    ##    browser()
    startX <- startDate - beforeDays
    endX <- endDate + moreDays    
    srcProj <- data.table(Date = seq.Date(startX, endDate + 52*7, by = 1))
    src <- merge(src, srcProj, by = 'Date', all=T)
    src[, reg := reg[1]]
    src$proj <- src[, rnbinom(length(fit), mu = fit, size  = nbsize)]
    if(is.null(xlim)) xlim <- c(startX,endX)
    if(doPlot) {
        src[Date < endDate & Date > startX, plot(Date, cases, type = 'h', bty = 'n', las = 2, xaxt = 'n', xlim=xlim, xlab = '', lwd = 3)]
        axis.Date(1, at = seq.Date(startX, endX, by = 7), labels = F)
        axis.Date(1, at = seq.Date(startX, endX, by = 14), format = '%b-%d', las = 2)    
        rect(startDate, 0, endDate, par('usr')[4], col = rgb(0,.5,0,.3), border=NA)
        title(main)
browser()
        src[Date > endDate, lines(Date, proj, type = 'h', lwd = 3, col = 'red')]
        src[Date > startDate, lines(Date, fit, lty = 1, col = 'dodger blue', lwd = 2)]
        rect(endDate, 0, endDate + moreDays, par('usr')[4], col = rgb(0.5,0,0,.3), border=NA)
        text(startDate+include_interval/2, par('usr')[4], 'fitting\n window', pos=1)
        text(endDate+moreDays/2, par('usr')[4], 'forecasting \nwindow', pos=1)
    }
    return(src)
})

createHazTraj <- function(fits, nbsize = .8, propInTrial = .03, numClus = 20, clusSize = 300, weeks = T) {
    hazTList <- NULL
    for(cc in 1:numClus) {
        fit <- fits[[sample(regs, 1)]]
        src <- forecast(fit, doPlot = F)
        src$day <- src[, as.numeric(Date - Sys.Date())]
        lastDataDay <- src[max(which(!is.na(src$cases))), day]
        src[day < lastDataDay & is.na(cases), cases := 0] ## fill in over interval without reporting
        src$haz <- src$cases
        src[day > lastDataDay, haz := proj]
        src[, haz := haz * propInTrial / clusSize]
        src[day > -60, list(day, haz)]
        if(weeks) {
            src$week <- src[, floor(day/7)]
            src <- src[, list(haz = mean(haz, na.rm=T)), week] ## get mean daily hazard by week
        }
        src$cluster <- cc
        hazTList[[cc]] <- as.data.frame(src)
    }
    hazT <- rbindlist(hazTList)
    setnames(hazT, 'haz', 'clusHaz')
    hazT$day <- hazT[, week*7]
    hazT <- arrange(hazT, day, cluster)
    return(hazT[, list(cluster, day, clusHaz)])
}


