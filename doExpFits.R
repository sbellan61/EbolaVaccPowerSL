## Steve Bellan 2015
if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
source('ExpFit.R'); require(data.table)
## truncDate <- as.Date('2015-02-09')
## sl <- sl[Date < truncDate]
xlim <- as.Date(c('2014-09-15','2015-05-01'))
 
regs <- levels(sl$reg)

## Fit to current incidence trends
include_interval <- 60 ## how many days back to include in exponential decay fit
minCases <- 30 ## if smaller than this then fit back to peak in subnational unit
fits <- NULL
for(rr in regs) fits[[rr]] <- doProj(sl[reg==rr], include_interval = include_interval, minCases = minCases, ll='exp_nbinom_ll', verbose=19)

nbsizeS <- unlist(lapply(fits, function(rr) rr$fit$par["nbsize"]))
mean(nbsizeS) ## 1.2

## Show fits and one simulated projection by subregion
##pdf('Figures/forecasted Paneled SL cleaned subnational data fromMax.pdf',  w = 10, h = 8)
jpeg('Figures/forecasted Paneled SL cleaned subnational data fromMax.jpg',  w = 10, h = 8, units='in',res=200)
nbsize <- 1.2 ## NULL
par(lwd=1, 'ps' = 12, mar = c(5,3,1.5,.5),mfrow = c(4,4))
regs <- sl[,unique(reg)]
srcs <- NULL
for(rr in regs) srcs[[rr]] <- forecast(fits[[rr]], main = rr, nbsize = nbsize, xlim = xlim,verbose = 19)
graphics.off()
srcProj <- rbindlist(srcs)

## Show simulated hazards from fits
pdf('Figures/example hazT.pdf')
for(jj in 1:10) {
    par(mar=c(5,5,2,.5), 'ps'=12, mgp = c(4,1,0))
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-month', main='', bty = 'n', las = 1,
         xlim = c(-5, 35), ylim = c(0,.02))
    ht <- createHazTraj(fits, numClus = 20, trialStartDate = as.Date('2015-02-01')) ## start date works, can test here
    ht[, lines(day/7,clusHaz*30, col = rainbow(20)[cluster], type = 'l', lwd = 2), cluster]
}
graphics.off()

save(fits, regs, file = 'data/createHT.Rdata')
