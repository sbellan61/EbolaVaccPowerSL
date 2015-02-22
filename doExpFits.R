## Steve Bellan 2015
if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
source('ExpFit.R'); require(data.table)

## allSL <- sl[,list(cases=sum(cases,na.rm=T)),Date]
## allSL <- allSL[,list(Date,cases)]

regs <- levels(sl$reg)

## Fit to current incidence trends
include_interval <- 60 ## how many days back to include in exponential decay fit
minCases <- 30 ## if smaller than this then fit back to peak in subnational unit
fits <- NULL
for(rr in regs) fits[[rr]] <- doProj(sl[reg==rr], include_interval = include_interval, minCases = minCases)

# sapply(fits, function(rr) rr$fit$par["nbsize"])

## Show fits and one simulated projection by subregion
pdf('Figures/forecasted Paneled SL cleaned subnational data fromMax.pdf',  w = 10, h = 8)
nbsize <- .8
par(lwd=1, 'ps' = 12, mar = c(5,3,1.5,.5),mfrow = c(4,4))
regs <- sl[,unique(reg)]
srcs <- NULL
xlim <- as.Date(c('2014-09-15','2015-05-01'))
for(rr in regs) srcs[[rr]] <- forecast(fits[[rr]], main = rr, nbsize = nbsize, xlim = xlim)
graphics.off()
srcProj <- rbindlist(srcs)

## Show simulated hazards from fits
pdf('Figures/example hazT.pdf')
for(jj in 1:10) {
    par(mar=c(5,5,2,.5), 'ps'=12, mgp = c(4,1,0))
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-month', main='', bty = 'n', las = 1,
         xlim = c(-5, 35), ylim = c(0,.02))
    ht <- createHazTraj(fits, numClus = 20)
    ht[, lines(day/7,clusHaz*30, col = rainbow(20)[cluster], type = 'l', lwd = 2), cluster]
}
graphics.off()

save(fits, regs, file = 'data/createHT.Rdata')

## Get distribution of decay rates
decayRates <- unlist(lapply(fits, function(x) as.numeric(x$fit$par['decay_rate'])))
wdr <- 1-(1-decayRates)^7

## Fit normal distribution function to decay rates
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

