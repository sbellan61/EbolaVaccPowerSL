if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)

pdf('example default cluster hazards.pdf')
hist(simTrial(makeParms('SWCT'))$pop[idByClus==1,clusHaz*30], col = 'black', breaks = 20, xlab = 'mean hazard by cluster (per month) ', 
     main ='distribution of mean hazard by cluster for 20 clusters')
graphics.off()

pdf('example individual hazard distr.pdf')
par(mfrow=c(5,4), mar = c(4,1,1,.5), oma = c(0,0,1.5,0))
tpop <- simTrial(makeParms('SWCT'))$pop
xlim <- c(0, max(tpop[,indivHaz]*30))
for(ii in 1:20) hist(tpop[cluster==ii,indivHaz*30], col = 'black', breaks = 20, xlab = 'individual hazard (per month)', main = '',
                     yaxt='n', xlim = xlim)
mtext('variation between clusters', 3,0, T)
graphics.off()

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
print('Done')
save(swctSims, rctSims, frctSims, crctSims, file = 'sims.Rdata')

pdf('stopping time results.pdf')
par(mfrow=c(4,1), mar = c(4,4,1,0))
xlim <- c(0,200)
for(typ in c('swct','rct','frct','crct')) {
    tempout <- as.data.table(get(paste0(typ,'Sims')))
    hist(tempout[,stopDay], breaks = seq(0, 10^4, by = 10), col = 'black', xlab ='stopping time', main = typ, xlim = xlim)
    abline(v=mean(tempout[,stopDay]), col = 'red')
}
graphics.off()
