if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)

nsims <- 100
swctSims <- simNtrials(parms=makeParms('SWCT'), N = nsims)
rctSims <- simNtrials(parms=makeParms('RCT'), N = nsims)
frctSims <- simNtrials(parms=makeParms('RCT', delayUnit = 3.5), N = nsims)
crctSims <- simNtrials(parms=makeParms('CRCT'), N = nsims)

## examine study designs
select(simTrial(makeParms('SWCT', numClus = 4, clusSize = 4))$pop, cluster, vaccDay)
select(simTrial(makeParms('RCT', numClus = 4, clusSize = 4))$pop, cluster, vaccDay)
select(simTrial(makeParms('RCT', numClus = 4, clusSize = 4, delayUnit = 3.5))$pop, cluster, vaccDay) ## RCT limited by
select(simTrial(makeParms('CRCT', numClus = 4, clusSize = 4))$pop, cluster, vaccDay)

par(mfrow=c(4,1), mar = c(4,4,1,0))
for(typ in c('swct','rct','frct','crct')) {
    tempout <- get(paste0(typ,'Sims'))
    hist(tempout[,stopDay], breaks = seq(0, max(tempout[,stopDay])+10, by = 10), col = 'black', xlab ='stopping time', main = typ)
    abline(v=mean(tempout[,stopDay]), col = 'red')
}
