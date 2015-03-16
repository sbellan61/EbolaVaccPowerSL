## Steve Bellan 2015
if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
source('ExpFit.R'); require(data.table)
## truncDate <- as.Date('2015-02-09')
## sl <- sl[Date < truncDate]
xlim <- as.Date(c('2014-09-15','2015-05-01'))
 
regs <- levels(sl$reg)

## Fit to current incidence trends
fits <- NULL
for(rr in regs) fits[[rr]] <- doProj(sl[reg==rr], ll='exp_nbinom_ll')

## Show fits and one simulated projection by subregion (4 districts only)
for(ii in 1:2) {
    if(ii==2) png('Figures/Fig 1 - comb.png',  w = 6.5, h = 3.5, units='in',res=200)
    if(ii==1) pdf('Figures/Fig 1 - comb.pdf',  w = 6.5, h = 3.5)
    nbsize <- 1.2 ## NULL
    par(lwd=1, 'ps' = 10, mar = c(1,3,2,.5),mfrow = c(2,3), oma = c(2,1.5,0,0))
    regs <- sl[,unique(reg)]
    srcs <- NULL
    regDo <- c('Kono','PortLoko', 'WesternAreaUrban', 'WesternAreaRural')
    labDo <- c('Kono','Port Loko', 'Western Area Urban', 'Western Area Rural')
    for(ri in 1:4) srcs[[regDo[ri]]] <- forecast(fits[[regDo[ri]]], main = labDo[ri], nbsize = nbsize
                                                 , xlim = xlim, xticks = ri>2)
    mtext('new cases', 2, 0, outer=T)
    ## Show simulated hazards from fits
    set.seed(8)
    #par(mar=c(3,1,2,.5), 'ps'=12, mgp = c(4,1,0), mfrow = c(1,2), oma = c(0,4,0,0))
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
    ##}
    title('(A)')##. Cluster-level hazards')
########## individual heterogneiety
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
    ## xlim <- c(.05, 20)
    ## xx <- exp(pretty(log(xlim), 1000))
    ## yy <- dnorm(log(xx), 0,  1)
    ## plot(xx,yy, type = 'h', xlab = 'risk factor', ylab='', main='', 
    ##      bty = 'n', xlim = xlim, log='x',ylim = c(0,.6), axes=F)
    ## ##labs <- c(expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0),expression(10^1), expression(10^2), expression(10^3))
    ## axis(1, at = c(.1,.5,1,2,10), las = 2)#, lab = c(expression(1/10),expression(1/2),1,2,10), las = 1)
    title('(B)') ##. Individual-level variation \naround cluster mean')
    mtext('hazard per person-month', 2, 3, outer = T, at = .6)
graphics.off()
}
