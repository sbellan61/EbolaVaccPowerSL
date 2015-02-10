## From https://github.com/ICI3D/ebola-forecast/ (Carl Pearson)
library(boot); library(ggplot2)
load(file='cleanSLData.Rdata')

param.xform <- function(x) c(x['rate_0'], decay_rate = as.numeric(inv.logit(x['logitdecay_rate'])))
param.xform(c(rate_0=30, logitdecay_rate=logit(0.01)))

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
exp_pois_ll <- llgenerator(exp_decay, pois_cases)

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

params <- function(src, censor_interval, include_interval, param.guess = c(rate_0=30, logitdecay_rate=logit(0.01)), browse = F) {
    if(browse) browser()
    end_date <- max(src$Date) - censor_interval
    start_date <- end_date - include_interval
    slice <- src[(start_date <= src$Date) & (src$Date <= end_date),] ## data.table doesn't work with shiny for some reason?
    slice$t <- as.numeric(slice$Date - min(slice$Date))
    if(sum(slice$int>1)>0) print('not always daily incidence')
    res <- optim(param.guess, exp_pois_ll, gr=NULL, cases = slice$cases, times = slice$t, control=list(fnscale=-1, parscale=c(1,1/10)), hessian=TRUE)
    vcmat <- solve(-res$hessian)
    res$par <- param.xform(res$par)
    cis <- with(res,{
        offset <- 1.96*sqrt(vcmat[2,2])
        par["decay_rate"] + c(low = offset, est=0, hi = -offset) # hi b = hi decay = low time
    })
    list(t0 = start_date, par=res$par, cis=cis)
}

allSL <- sl[,list(cases=sum(cases,na.rm=T)),Date]
allSL <- allSL[,list(Date,cases)]

doProj <- function(src, forecast_date = tail(src$Date,1), censor_interval = 0, include_interval = 30, minCases = 15, 
                      beforeDays = 30, moreDays = 90, minDecayRate = .01) {
    endDate <- max(src$Date) - censor_interval
    startDate <- endDate - include_interval
    startX <- startDate - beforeDays
    ## minimum cases in last interval to fit exponential, otherwise increase to max # cases
    fromMax <- ifelse(src[Date > max(Date)-include_interval, sum(cases)] < minCases, T, F)
    if(fromMax) startDate <- src[which.max(cases), Date]
    include_interval <- endDate-startDate
    fit_ref <- params(src, censor_interval, include_interval)
    fit_ref$par['decay_rate'] <- max(fit_ref$par['decay_rate'], minDecayRate)
    src$fit <- expcurve(fit_ref, date_zero = as.numeric(startDate))(as.numeric(src$Date))
    return(list(fit = fit_ref, src = src[,list(Date,reg, cases,fit)], startDate = startDate, endDate = endDate,
                moreDays=moreDays, beforeDays=beforeDays, include_interval=include_interval))
}

forecast <- function(fit, main=NULL, nbsize = .9, doPlot = T) with(fit, {
    ##    browser()
    startX <- startDate - beforeDays
    endX <- endDate + moreDays    
    srcProj <- data.table(Date = seq.Date(startX, startX + 365, by = 1))
    src <- merge(src, srcProj, by = 'Date', all=T)
    src[, reg := reg[1]]
    src$fit <- expcurve(fit, date_zero = as.numeric(startDate))(as.numeric(src$Date))
    src$proj <- src[, rnbinom(length(fit), mu = fit, size  = nbsize)]
    if(doPlot) {
        src[Date < endDate & Date > startX, plot(Date, cases, type = 'h', bty = 'n', las = 2, xaxt = 'n', xlim=c(startX,endX), xlab = '', lwd = 3)]
        axis.Date(1, at = seq.Date(startX, endX, by = 7), labels = F)
        axis.Date(1, at = seq.Date(startX, endX, by = 14), format = '%b-%d', las = 2)    
        rect(startDate, 0, endDate, par('usr')[4], col = rgb(0,.5,0,.3), border=NA)
        title(main)
        src[Date > endDate, lines(Date, proj, type = 'h', lwd = 3, col = 'red')]
        src[, lines(Date, fit, lty = 1, col = 'dodger blue', lwd = 2)]
        rect(endDate, 0, endDate + moreDays, par('usr')[4], col = rgb(0.5,0,0,.3), border=NA)
        text(startDate+include_interval/2, par('usr')[4], 'fitting\n window', pos=1)
        text(endDate+moreDays/2, par('usr')[4], 'forecasting \nwindow', pos=1)
    }
    return(src)
})

createHazTraj <- function(fit, nbsize = .8, propInTrial = .03, clusSize = 300, weeks = T) {
    src <- forecast(fit, doPlot = F)
    src$day <- src[, as.numeric(Date - Sys.Date())]
    src$haz <- src$cases
    src[day > -20, haz := proj]
    src[, haz := haz * propInTrial / clusSize]
    src[day > -60, list(day, haz)]
    if(weeks) {
        src$week <- src[, floor(day/7)]
        src <- src[, list(haz = mean(haz, na.rm=T)), week] ## get mean daily hazard by week
    }
    return(src)
}
 
pdf('Figures/example hazT.pdf')
for(jj in 1:10) {
    par(mar=c(5,5,2,.5), 'ps'=12, mgp = c(4,1,0))
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-month', main='', bty = 'n', las = 1,
         xlim = c(0, 35), ylim = c(0,.015))
    rrs <- factor(NULL, levels = regs)
    for(ii in 1:20) {
        rr <- sample(regs, 1)
        rrs[ii] <- rr
        ht <- createHazTraj(fits[[rr]])
        lines(ht[,list(week,haz*30)], type = 'l', col = rainbow(20)[ii], lwd = 2)
    }
    legend('topright', leg = rrs, col = rainbow(20), lwd = 1, ncol = 3, cex = .8, bty = 'n')
}
graphics.off()

forecastInc(fits)

plotFit(fits[[1]])

## By subregion

fits <- NULL
for(rr in regs) fits[[rr]] <- doProj(sl[reg==rr], include_interval = 60, minCases = 30)

pdf('../Figures/forecasted Paneled SL cleaned subnational data fromMax.pdf',  w = 10, h = 8)
nbsize <- .8
par(lwd=1, 'ps' = 12, mar = c(5,3,1.5,.5),mfrow = c(4,4))
regs <- sl[,unique(reg)]
srcs <- NULL
for(rr in regs) srcs[[rr]] <- forecast(fits[[rr]], main = rr, nbsize = nbsize)
graphics.off()
srcProj <- rbindlist(srcs)

## Get distribution of decay rates
decayRates <- unlist(lapply(fits, function(x) as.numeric(x$fit$par['decay_rate'])))
wdr <- 1-(1-decayRates)^7

## Fit normal distribution function
wdrFit <- fitdistr(logit(wdr), dnorm, start = list(mean = -1, sd = .5))$estimate

pdf('../Figures/histogram of logit weekly decay rates.pdf', w = 5, h = 5.5)
par('ps'=12)
hist(logit(wdr), breaks = 10, col = 'black', xlab = 'logit(weekly decay rate)', main = 'variation in decay rate', freq = F)
xs <- seq(-3,3, by = .1)
ys <- dnorm(xs, wdrFit['mean'], wdrFit['sd'])
lines(xs, ys, lwd = 2, col = 'dodger blue')
graphics.off()

pdf('../Figures/histogram of weekly decay rates.pdf', w = 5, h = 5.5)
par('ps'=12)
h1 <- hist(wdr, breaks = 10, col = 'black', xlab = 'weekly decay rate', main = 'variation in decay rate', xlim = c(0, .7), freq = F)
xs <- seq(-3,3, by = .1)
ys <- dnorm(xs, wdrFit['mean'], wdrFit['sd'])
ys <- ys*max(h1$density)/max(ys)
lines(inv.logit(xs), ys, lwd = 2, col = 'dodger blue')
graphics.off()

decayRateFXN <- function(n, wdrparms = wdrFit, unit = 'week') {
    logitwdr <- rnorm(n, wdrparms['mean'], wdrparms['sd'])
    rand <- inv.logit(logitwdr)
    if(unit == 'day') {
        rand <- (1 - rand)^(1/7)
        rand <- 1-rand
    }
    return(rand)
}

decayRateFXN(10, unit='day')
save(wdrFit, decayRateFXN, file = 'wdrRNG.Rdata')

