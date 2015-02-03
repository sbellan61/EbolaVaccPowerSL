outdir <- 'Figures'

plotHazT <- function(parms, flnm = NULL, browse=F, main='', ...) with(parms, {
    if(!is.null(flnm)) pdf(file.path(outdir, paste0(flnm,'.pdf')), ...)
    plot(0,0, type = 'n', xlab = 'weeks', ylab = 'hazard per person-year', main=main, bty = 'n', las = 1,
         xlim = c(0, (maxInfectDay+hazIntUnit)/7), ylim = c(0,max(popH$clusHaz)/yearToDays))
    popH[idByClus==1, lines(day/7, clusHaz/yearToDays, type = 's', col = rainbow(numClus)[cluster]), by = cluster]
    if(!is.null(flnm)) graphics.off()
})

p1 <- simTrial(makeParms('RCT',small=F, ord='TU', delayUnit = 0, clusSize=1, numClus = 20))
nm <- 'default variation in EVD hazard \nfor trial participants'
plotHazT(p1, flnm='hazard trajectories' , main=nm )


plotTrialRollout <- function(parms, flnm = NULL, browse=F, main='', ...) with(parms, {
})
