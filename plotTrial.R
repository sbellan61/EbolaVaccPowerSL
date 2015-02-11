if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
outdir <- 'Figures'

plotHazT <- function(parms, flnm = NULL, browse=F, main='', ymax = NULL, ...) with(parms, {
    if(!is.null(flnm)) pdf(file.path(outdir, paste0(flnm,'.pdf')), ...)
    par(mgp=c(4,1,0), mar = c(5,5,2,.5))
    if(is.null(ymax)) ymax <- popH[, max(clusHaz)] * 30
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-month', main=main, bty = 'n', las = 1,
         xlim = c(0, (maxInfectDay+hazIntUnit)/7), ylim = c(0,ymax)) #max(popH$clusHaz)*30))
    popH[idByClus==1, lines(day/7, clusHaz*30, type = 'l', col = rainbow(numClus)[cluster]), by = cluster]
    if(!is.null(flnm)) graphics.off()
})

p1 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, propInTrial = .03, numClus = 20), br = F)
p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazSL = F, weeklyDecay = .9, cvWeeklyDecay = .5, cvClus = 1.5, numClus = 20), br = F)

plotHazT(p1, flnm='SL hazard trajectories' , main='mean cluster hazards from SL forecasts\n (3% district cases in trial)', ymax = .01)
plotHazT(p2, flnm='sim hazard trajectories' , main='mean cluster hazards from\n (weekly decay = .9, decayCV = .5, baselineCV=1.5)', ymax = .01)


plotTrialRollout <- function(parms, flnm = NULL, browse=F, main='', ...) with(parms, {
})
