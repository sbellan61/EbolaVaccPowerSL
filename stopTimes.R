if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)

## examine study designs
select(simTrial(makeParms('SWCT', numClus = 4, clusSize = 4))$pop, cluster, vaccDay)
select(simTrial(makeParms('RCT', numClus = 4, clusSize = 4))$pop, cluster, vaccDay)
select(simTrial(makeParms('RCT', numClus = 4, clusSize = 4, delayUnit = 3.5))$pop, cluster, vaccDay) ## RCT limited by
select(simTrial(makeParms('CRCT', numClus = 4, clusSize = 4))$pop, cluster, vaccDay)

nsims <- 500
system.time(swctSims <- simNtrials(parms=makeParms('SWCT'), N = nsims, verbose=0))
rctSims <- simNtrials(parms=makeParms('RCT'), N = nsims)
frctSims <- simNtrials(parms=makeParms('RCT', delayUnit = 3.5), N = nsims)
crctSims <- simNtrials(parms=makeParms('CRCT'), N = nsims)
save(swctSims, rctSims, frctSims, crctSims, file = 'sims.Rdata')

pdf('stopping time results.pdf')
par(mfrow=c(4,1), mar = c(4,4,1,0))
xlim <- c(0,200)
for(typ in c('swct','rct','frct','crct')) {
    tempout <- get(paste0(typ,'Sims'))
    hist(tempout[,stopDay], breaks = seq(0, 10^4, by = 10), col = 'black', xlab ='stopping time', main = typ, xlim = xlim)
    abline(v=mean(tempout[,stopDay]), col = 'red')
}
graphics.off()
