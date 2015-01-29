if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)

pdf('distribution of mean cluster hazards.pdf')
parms <- makeParms()
samp <- reParmRgamma(10^5, parms$mu, parms$varClus)
breaks <- 100
hist(samp/monthToDays, col = 'black', xlab = 'hazard / month', main = 'distribution of cluster hazards', 
     xlim = c(0, .05), breaks=breaks)
title(main=paste('average hazard = ', signif(parms$mu/monthToDays,2)), line = -2)
graphics.off()

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
simTrial(makeParms('SWCT', numClus = 4, clusSize = 4))$pop[ , list( cluster, vaccDay)]
simTrial(makeParms('RCT', numClus = 4, clusSize = 4))$pop[ , list( cluster, vaccDay)]
simTrial(makeParms('FRCT', numClus = 4, clusSize = 4))$pop[ , list( cluster, vaccDay)]
simTrial(makeParms('CRCT', numClus = 4, clusSize = 4))$pop[ , list( cluster, vaccDay)]
