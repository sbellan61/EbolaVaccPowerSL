####################################################################################################
## Perform & plot Sierra Leone district-level exponential model fits along with incidence
## projections. Plot example cluster-level hazard trajectories and individual-level variation around
## cluster means.
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

source('ExpFit.R'); require(data.table)
xlim <- as.Date(c('2014-09-15','2015-05-01'))
regs <- levels(sl$reg)

## Fit to current incidence trends
fits <- NULL
for(rr in regs) fits[[rr]] <- doProj(sl[reg==rr], ll='exp_nbinom_ll')

## Overdispersion parameters
nbsizeS <- unlist(lapply(fits, function(rr) rr$fit$par["nbsize"]))
mean(nbsizeS) ## 1.2

####################################################################################################
## Figure 1
####################################################################################################
## Show fits and one simulated projection by subregion (4 districts only)
for(ii in 1:2) {
    if(ii==2) png('Figures/Figure 1 - forecasted district data.png',  w = 6.5, h = 3.5, units='in',res=200)
    if(ii==1) pdf('Figures/Figure 1 - forecasted district data.pdf',  w = 6.5, h = 3.5)
    set.seed(2)
    nbsize <- 1.2 ## NULL
    par(lwd=1, 'ps' = 10, mar = c(1,3,2,.5),mfrow = c(2,2), oma = c(2,1.5,0,0))
    regs <- sl[,unique(reg)]
    srcs <- NULL
    regDo <- c('Kono','PortLoko', 'WesternAreaUrban', 'WesternAreaRural')
    labDo <- c('Kono','Port Loko', 'Western Area Urban', 'Western Area Rural')
    for(ri in 1:4) srcs[[regDo[ri]]] <- forecast(fits[[regDo[ri]]], main = labDo[ri], nbsize = nbsize
                                               , xlim = xlim, xticks = ri>2)
    mtext('new cases', 2, 0, outer=T)
    graphics.off()
}

####################################################################################################
## Figure S11
####################################################################################################
## Show fits and one simulated projection by subregion
for(ii in 1:2) {
    set.seed(1)
    if(ii==1) jpeg('Figures/Figure S11 - forecasted Paneled SL cleaned subnational data fromMax.jpg',  w = 10, h = 8, units='in',res=200)
    if(ii==2) pdf('Figures/Figure S11 - forecasted Paneled SL cleaned subnational data fromMax.pdf',  w = 10, h = 8)
    nbsize <- 1.2 ## NULL
    par(lwd=1, 'ps' = 10, mar = c(5,3,1.5,.5),mfrow = c(4,4))
    regs <- sl[,unique(reg)]
    srcs <- NULL
    for(ri in 1:length(regs)) {
        rr <- regs[ri]
        srcs[[rr]] <- forecast(fits[[rr]], main = rr, nbsize = nbsize, xlim = xlim)#, xticks = ri>12)
    }
    graphics.off()
    srcProj <- rbindlist(srcs)
}
 
####################################################################################################
## Figure 2
####################################################################################################
## Show simulated hazards from fits
for(ii in 1:2) {
    if(ii==2) png('Figures/Figure 2 - example hazT.png', w = 6.5, h = 4, units='in', res = 200)
    if(ii==1) pdf('Figures/Figure 2 - example hazT.pdf', w = 6.5, h = 4)
    set.seed(8)
    par(mar=c(3,1,2,.5), 'ps'=12, mgp = c(4,1,0), mfrow = c(1,2), oma = c(0,4,0,0))
    ## Cluster-level heterogeneity
    ht <- createHazTrajFromSLProjection(fits, numClus = 20, trialStartDate = as.Date('2015-02-01')) ## start date works, can test here
    xlim <- ht[day>=-30 & day < 150, range(Date)]
    plot(xlim,c(0,0), type = 'n', xlab = '', ylab = 'hazard per person-month', main='', bty = 'n', las = 2,
         xlim = as.Date(xlim), ylim = c(0,.03), axes = F)
    yticks <- seq(0,.03,by=.005)
    axis(2, yticks, las = 2)
    mnth <- seq.Date(as.Date('2015-01-01'),as.Date('2015-07-01'), by ='month')
    mnth.mids <- seq.Date(as.Date('2015-01-15'),as.Date('2015-07-15'), by ='month')
    axis.Date(1, at = mnth,  format = '%b-%d', label = NA)
    axis.Date(1, at = mnth.mids,  format = '%b', las = 2, tck=0)
    ht[, lines(Date,clusHaz*30, col = rainbow(20)[cluster], type = 'l', lwd = 2), cluster]
    title('(A)') 
    ## Individual heterogneiety
    plot(xlim,c(0,0), type = 'n', xlab = '', ylab = '', main='', bty = 'n', las = 2,
         xlim = as.Date(xlim), ylim = c(0,.03), axes=F)
    axis(2, yticks, lab = NA)
    axis.Date(1, at = mnth,  format = '%b-%d', label = NA)
    axis.Date(1, at = mnth.mids,  format = '%b', las = 2, tck=0)
    nrand <- 2
    clsh <- 2
    rand <- data.table(indiv = 1:nrand, rf = qlnorm(c(.25,.75)), cluster = clsh)
    htshow <- ht[cluster==clsh]
    htshow <- merge(htshow, rand, by = 'cluster', allow.cartesian=T)
    htshow[, indivHaz := clusHaz* rf]
    coltr <- makeTransparent(rainbow(20)[clsh])
    polygon(htshow[indiv==1,c(Date,rev(Date))], c(htshow[indiv==1, indivHaz*30], rev(htshow[indiv==2, indivHaz*30])), 
            col = coltr , border = NA)
    ht[cluster==clsh, lines(Date,clusHaz*30, col = rainbow(20)[cluster], type = 'l', lwd = 2)]
    title('(B)') 
    mtext('hazard per person-month', 2, 3, outer = T, at = .6)
    graphics.off()
}
save(fits, regs, ht, file = 'data/createHT.Rdata')

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
