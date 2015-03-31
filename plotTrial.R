####################################################################################################
## Plot simulated hazard trajectories.
####################################################################################################
## Code base accompanying:
## 
## Bellan, SE, JRC Pulliam, CAB Pearson, DChampredon, SJ Fox, L Skrip, AP Galvani, M Gambhir, BA
## Lopman, TC Porco, LA Meyers, J Dushoff (2015). The statistical power and validity of Ebola
## vaccine trials in Sierra Leone: A simulation study of trial design and analysis. _Lancet
## Infectious Diseases_.
##
## Steve Bellan, March 2015
## License at bottom.
####################################################################################################
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','ExpFit.R'), source)
outdir <- 'Figures'

## Plot simulated hazard trajectories
plotHazT <- function(parms, flnm = NULL, browse=F, main='', ymax = NULL, ytickLab = T
                   , xlab = 'weeks', ylab = 'hazard per person-month', ...) 
    with(parms, {
        if(!is.null(flnm)) {
            pdf(file.path(outdir, paste0(flnm,'.pdf')), ...)
            par(mgp=c(4,1,0), mar = c(5,5,2,.5))
        }
        if(is.null(ymax)) ymax <- popH[, max(clusHaz)] * 30
        plot(0,0, type = 'n', xlab = xlab, ylab = ylab, main=main, bty = 'n', las = 1, axes = F,
             xlim = c(0, (maxInfectDay+hazIntUnit)/7), ylim = c(0,ymax)) #max(popH$clusHaz)*30))
        if(ytickLab) axis(2, seq(0, ymax, l =5), las = 1) else axis(2, seq(0, ymax, l =5), labels = NA)
        popH[idByClus==1, lines(day/7, clusHaz*30, type = 'l', col = rainbow(numClus)[cluster]), by = cluster]
        axis(1, at = c(0,12,24), las = 1)
        axis(1, at = c(6, 18), labels = NA)
        if(!is.null(flnm)) graphics.off()
    })

p1 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'SL', propInTrial = .03, numClus = 20))
p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'Phenom', weeklyDecay = .9, cvWeeklyDecay = .5, cvClus = 1.5, cvClusTime = 0.5, numClus = 20))

jpeg('Figures/sim Sl hazards.jpg', w = 4.5, h = 4, units='in', res = 200)
par('ps'=10, mgp = c(3,1,0), mar = c(5,4,2,.5))
plotHazT(p1, main='simulated hazards based \non Sierra Leone forecast', ymax = .02)
graphics.off()

pdf('Figures/sim hazards.pdf', w = 6.5, h = 4) ## show 10 simulations to see variation
for(ii in 1:10) {
    p1 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'SL', propInTrial = .03, numClus = 20))
    p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'Phenom', weeklyDecay = .9, cvWeeklyDecay = .5, cvClus = 1.8, cvClusTime = 0.5, numClus = 20))
    par(mfrow=c(1,2), 'ps'=10, mgp = c(3.5,1,0), mar = c(5,5,2,.5))
    plotHazT(p1, main='(A) simulated hazards based \non Sierra Leone forecasts', ymax = .02)
    plotHazT(p2, main='(B) simulated exponential \ndecay hazards', ymax = .02)
}
graphics.off()

####################################################################################################
## Figure S6
####################################################################################################
pdf('Figures/Figure S6 - phenom haz.pdf', w = 6.5, h = 5)
cvwd <- c(0,.5,1)
cvct <- c(0,.25,.5,1)
seed <- 2
par(mfrow=c(3,4), 'ps'=10, mgp = c(3.5,1,0), mar = c(2,1.5,2,.5), oma = c(1,6,1,0))
for(kk in 1:10) {
for(jj in 1:3) {
    for(ii in 1:4) {
        p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'Phenom', weeklyDecay = .9, 
                                 cvWeeklyDecay = cvwd[jj], cvClus = 1.8, cvClusTime = cvct[ii], numClus = 20, dontReordForPlot=T), seed = kk)
        plotHazT(p2, main='', ymax = .02, xlab = '', ylab = '', ytickLab = ii==1)
    }}
mtext('increasing variation around trend', 3, 0, T)
mtext('increasing variation in trend', 2, 4.5, T)
mtext('hazard per person-month', 2, 2, T)
mtext('weeks', 1, 0, T)
}
graphics.off()

####################################################################################################
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/sbellan61/EbolaVaccPowerSL
####################################################################################################
