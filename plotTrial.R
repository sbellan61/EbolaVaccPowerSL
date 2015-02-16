if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R','ExpFit.R'), source)
outdir <- 'Figures'

plotHazT <- function(parms, flnm = NULL, browse=F, main='', ymax = NULL, ...) with(parms, {
    if(!is.null(flnm)) {
        pdf(file.path(outdir, paste0(flnm,'.pdf')), ...)
        par(mgp=c(4,1,0), mar = c(5,5,2,.5))
    }
    if(is.null(ymax)) ymax <- popH[, max(clusHaz)] * 30
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-month', main=main, bty = 'n', las = 1,
         xlim = c(0, (maxInfectDay+hazIntUnit)/7), ylim = c(0,ymax)) #max(popH$clusHaz)*30))
    popH[idByClus==1, lines(day/7, clusHaz*30, type = 'l', col = rainbow(numClus)[cluster]), by = cluster]
    if(!is.null(flnm)) graphics.off()
})

p1 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, propInTrial = .03, numClus = 20), br = F)

jpeg('Figures/sim Sl hazards.jpg', w = 4.5, h = 4, units='in', res = 200)
par('ps'=10, mgp = c(3,1,0), mar = c(5,4,2,.5))
plotHazT(p1, main='simulated hazards based \non Sierra Leone forecast', ymax = .01)
graphics.off()

p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazSL = F, weeklyDecay = .9, cvWeeklyDecay = 1, cvClus = 1.5, numClus = 20), br = F)
plotHazT(p2, flnm='sim hazard trajectories' , main='mean cluster hazards from\n (weekly decay = .9, decayCV = .5, baselineCV=1.5)', ymax = .01)


jpeg('Figures/sim hazards.jpg', w = 6.5, h = 4, units='in', res = 200)
par(mfrow=c(1,2), 'ps'=10, mgp = c(3.5,1,0), mar = c(5,5,2,.5))
plotHazT(p1, main='(A) simulated hazards based \non Sierra Leone forecasts', ymax = .01)
plotHazT(p2, main='(B) simulated exponential \ndecay hazards', ymax = .01)
dev.off()

plotTrialRollout <- function(parms, flnm = NULL, browse=F, main='', ...) with(parms, {
})


showSeqStop <- function(resfull, flnm= NULL, ...) {
    with(resfull, {
        if(!is.null(flnm)) pdf(file.path('Figures', paste0(flnm,'.pdf')), ...)
        par(lwd=1.5, mar = c(3,6,1,5), mgp = c(4,1,0), mfrow = c(2,2), oma = c(2,0,0,0))
        plot(0,0, type = 'n', xlim = range(popH$day)/7, ylim = c(0, mu*7), las= 1, bty = 'n',
             ylab = 'weekly hazard', xlab = 'day')
        lines(weeklyAns[,list(stopDay/7, hazCXimmGrpEnd)], col = 'red')
        lines(weeklyAns[,list(stopDay/7, hazVXimmGrpEnd)], col = 'black')
        legend('topleft', c('vacc','cont','total','P value'), col = c('black', 'red','dark green','purple'),  bty = 'n', lwd=2,bg='white')
        ## person-year of observation
        plot(0,0, type = 'n', xlim = range(popH$day)/7, ylim = c(0, yearToDays*max(weeklyAns[, caseCXimmGrpEnd/hazCXimmGrpEnd],na.rm=T)), las= 1, bty = 'n',
             ylab = 'person-years', xlab = 'day')
        lines(weeklyAns[,list(stopDay/7, yearToDays*caseCXimmGrpEnd/hazCXimmGrpEnd)], col = 'red')
        lines(weeklyAns[,list(stopDay/7, yearToDays*caseVXimmGrpEnd/hazVXimmGrpEnd)], col = 'black', lty = 2)
        legend('topleft', c('vacc','cont','P value'), col = c('black', 'red','purple'),  bty = 'n', lwd=2,bg='white')
        ## number of cases
        plot(0,0, type = 'n', xlim = range(popH$day)/7, ylim = c(0, max(weeklyAns[, list(caseVXimmGrpEnd,caseCXimmGrpEnd)])), las= 1, bty = 'n',
             ylab = 'cases', xlab = 'day')
        lines(weeklyAns[,list(stopDay/7, caseCXimmGrpEnd)], col = 'red')
        lines(weeklyAns[,list(stopDay/7, caseVXimmGrpEnd)], col = 'black')
        lines(weeklyAns[,list(stopDay/7, caseVXimmGrpEnd+caseCXimmGrpEnd)], col = 'dark green')
        par(new=T)
        plot(weeklyAns[,list(stopDay, p)], col = 'purple', lty = 1, axes = F, ylab='', xlab='', type='l', ylim = c(0,1))
        axis(4, at = seq(0, 1, by = .05), lab = NA)
        axis(4, at = seq(0, 1, by = .1), las = 1)
        abline(h=.05, lty = 2)
        mtext('p value', 4, 3)
        ## vaccine efficacy estimate
        plot(0,0, type = 'n', xlim = range(popH$day)/7,
             ylim = c(-1, 1), las= 1, bty = 'n', ylab = 'vaccine efficacy', xlab = 'day')
        nmisg <- !weeklyAns[, is.na(lci) | is.na(uci)]
        polygon(c(weeklyAns$stopDay[nmisg], rev(weeklyAns$stopDay[nmisg]))/7,
                c(weeklyAns$lci[nmisg], rev(weeklyAns$uci[nmisg])), col = 'gray', border = NA)
        lines(weeklyAns[,list(stopDay/7, mean)], col = 'black')
        abline(h=0, lty = 2)
        mtext('week of trial', 1, 1, T)
        if(!is.null(flnm)) graphics.off()
    })
}
