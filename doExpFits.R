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

fits <- NULL
for(rr in regs) fits[[rr]] <- doProj(sl[reg==rr], ll='exp_nbinom_ll', verbose=19)

nbsizeS <- unlist(lapply(fits, function(rr) rr$fit$par["nbsize"]))
mean(nbsizeS) ## 1.2

## Show fits and one simulated projection by subregion (4 districts only)
png('Figures/Fig 1 - forecasted district data.png',  w = 6.5, h = 3.5, units='in',res=200)
nbsize <- 1.2 ## NULL
par(lwd=1, 'ps' = 10, mar = c(1,3,2,.5),mfrow = c(2,2), oma = c(3,1.5,0,0))
regs <- sl[,unique(reg)]
srcs <- NULL
regDo <- c('Kono','PortLoko', 'WesternAreaUrban', 'WesternAreaRural')
labDo <- c('Kono','Port Loko', 'Western Area Urban', 'Western Area Rural')
for(ri in 1:4) srcs[[regDo[ri]]] <- forecast(fits[[regDo[ri]]], main = labDo[ri], nbsize = nbsize
##                                             ,ylim = c(0,50)
                                             , xlim = xlim, xticks = ri>2, verbose = 19)
mtext('new cases', 2, 0, outer=T)
graphics.off()

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
