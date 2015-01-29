if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)

nsims <- 2

## NperCore <- 1
## system.time(
##     swctSims <- simNwrp(parms=makeParms('SWCT'), NperCore = 2, ncores = 1)
##     ## rctSims <- simNwrp(parms=makeParms('RCT'), NperCore = nsims)
##     ## frctSims <- simNwrp(parms=makeParms('RCT', delayUnit = 3.5), NperCore = nsims)
##     ## crctSims <- simNwrp(parms=makeParms('CRCT'), NperCore = nsims)
##     )

## makeDT <- function(x) data.table(1:10,21:30)

## mclapply(1:2, makeDT, mc.cores = 12)

## out <- mclapply(1:2, simNtrials, N = 1, mc.cores = 12)

system.time(swctSims <- simNtrials(parms=makeParms('SWCT'), N = nsims, verbose=0))
## rctSims <- simNtrials(parms=makeParms('RCT'), N = nsims)
## frctSims <- simNtrials(parms=makeParms('RCT', delayUnit = 3.5), N = nsims)
## crctSims <- simNtrials(parms=makeParms('CRCT'), N = nsims)


## save(swctSims, rctSims, frctSims, crctSims, file = 'sims.Rdata')

## pdf('stopping time results.pdf')
## par(mfrow=c(4,1), mar = c(4,4,1,0))
## xlim <- c(0,200)
## for(typ in c('swct','rct','frct','crct')) {
##     tempout <- as.data.table(get(paste0(typ,'Sims')))
##     hist(tempout[,stopDay], breaks = seq(0, 10^4, by = 10), col = 'black', xlab ='stopping time', main = typ, xlim = xlim)
##     abline(v=mean(tempout[,stopDay]), col = 'red')
## }
## graphics.off()
