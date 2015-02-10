if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
outdir <- 'Figures'

plotHazT <- function(parms, flnm = NULL, browse=F, main='', ...) with(parms, {
    if(!is.null(flnm)) pdf(file.path(outdir, paste0(flnm,'.pdf')), ...)
    par(mgp=c(4,1,0), mar = c(5,5,2,.5))
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-month', main=main, bty = 'n', las = 1,
         xlim = c(0, (maxInfectDay+hazIntUnit)/7), ylim = c(0,.015)) #max(popH$clusHaz)*30))
    popH[idByClus==1, lines(day/7, clusHaz*30, type = 's', col = rainbow(numClus)[cluster]), by = cluster]
    if(!is.null(flnm)) graphics.off()
})

p1 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=1, numClus = 20), br = F)
nm <- 'default variation in EVD hazard \nfor trial participants'
plotHazT(p1, flnm='hazard trajectories' , main=nm )


plotTrialRollout <- function(parms, flnm = NULL, browse=F, main='', ...) with(parms, {
})
