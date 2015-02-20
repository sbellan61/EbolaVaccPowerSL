if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)
##load(file=file.path('BigResults','powFin.Rdata'))

thing <- 'SLSims5'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))

pf[delayUnit==0, ord:='simultaneous instant']
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.03 & mod=='coxME']
setnames(pf,c('ord','trial'),c('order','design'))
levels(pf$order)[2] <- 'time-updated'
pf[, pit:=factor(paste0(propInTrial*100,'%'), levels = c('3%','5%','10%'), ordered=T)]

baseMods <- c('Cox PH Frailty'
              , 'Poisson GLM\n no cluster effects'
              , 'Poisson GLM \nwith fixed effects by cluster')

pf$modLb <- pf$mod
levels(pf$modLb) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))

##################################################
## Power
labs <- c('','log')
for(jj in 1:2) {
    pdf(paste0('Figures/',labs[jj], 'power SL ALL.pdf'), w = 8, h = 5) ##, units = 'in', res = 200)
    for(modTmp in levels(pf$modLb)) {
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[modLb==modTmp], aes(vaccEff, vaccGoodNAR, colour=design, linetype=order)) + 
            scale_x_continuous(labels = formatC, limits=c(0,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
                theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
                      axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
                      panel.margin = unit(2, "lines")) +
                          xlab('vaccine efficacy') + ylab('probability of rejecting null') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                  geom_line(size=1) + facet_wrap(~propInTrial) + 
                                      ggtitle('power by expected % of district-level cases in trial population') +
                                          geom_hline(yintercept=.025, color='black', linetype='dotted')
        if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), minor_breaks = seq(0,1,.05))
        if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, limits=c(0.01,.9), breaks = c(.01,.025,.05,.1,.2,.5,.8))
        print(p.tmp)
        grid.text(modTmp,x=unit(.85,"npc"),y=unit(0.85,"npc"), gp=gpar(fontsize=10))
    }
    graphics.off()
}

##################################################
## Type I
labs <- c('','log')
jj <- 1
for(jj in 1:2) {
pdf(paste0('Figures/',labs[jj],'Type I SL.pdf'), w = 8, h = 8) ##, units = 'in', res = 200)
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[design %in% c('SWCT','FRCT') & vaccEff==0], aes(propInTrial, stoppedNAR, colour=design, linetype=order)) + 
    scale_x_continuous(labels = formatC, limits=c(.03,.1), breaks = c(.03,.05,.1)) +  
##    scale_y_continuous(labels = formatC, limits=c(0,.3)) + #, minor_breaks = seq(0,1,.05)) +
    theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
          panel.margin = unit(2, "lines")) +
    xlab('% of distrinct-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~modLb) + 
geom_hline(yintercept=.05, color='black', linetype='dotted', size = 1.2)
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC)#, limits=c(0,.3))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .2, .3, .4))#,  limits=c(.01,.3))
print(p.tmp)
graphics.off()
}

##################################################
## Cases
pdf(paste0('Figures/Cases SL ALL.pdf'), w = 8, h = 5) ##, units = 'in', res = 200)
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[modLb==modTmp], aes(vaccEff, caseTot, colour=design, linetype=order)) + 
    scale_x_continuous(labels = formatC, limits=c(0,1), breaks = pf[,unique(vaccEff)],minor_breaks=NULL) +  
    theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5)) +
    xlab('vaccine efficacy') + ylab('# EVD cases') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit) + 
    ggtitle('# of cases in trial by expected % of district-level cases in trial population') +
    geom_hline(yintercept=.025, color='black', linetype='dotted')
p.tmp <- p.tmp + scale_y_continuous(labels = formatC) #
print(p.tmp)
graphics.off()


####################################################################################################
## numsims
pdf(paste0('Figures/NSIMS SL fp.pdf'), w = 8, h = 6) ##, units = 'in', res = 200)
par(lwd=2, mar = c(5,5,3,.5), oma = c(1,1,1,0))
layout(matrix(1:2,1,2), w = c(3,1))
subs <- pf[, mod==modtypes[1] & !(delayUnit==0 & ord!='none') & 
               ((trial=='SWCT' & ord=='none') | trial %in% c('FRCT','RCT','CRCT')) 
               & vaccEff ==0 ]
powTmp <- pf[subs]
powTmp[,trial:=factor(trial)]
plot(0,0, type = 'n', xlab = '', ylab = '', xlim = c(0.03,.1), ylim = c(0,1200), bty = 'n', main = '', las = 1, xaxt='n')
axis(1, at = c(0,.03, .05, .1))
powTmp[,
       {
           lty <- which(c('none','TU')==ord)
           if(delayUnit==0) lty <- 3
           lines(propInTrial, nsim, col = cols[tri==trial, col], lty = lty, type = 'b')
       },
       by = list(trial, ord, delayUnit > 0)]
par(mar=rep(0,4))
plot.new()
legend('topleft', leg=cols[tri %in% powTmp[,levels(trial)],tri], col = cols[tri %in% powTmp[,levels(trial)],col], lwd = 2, bty = 'n')
legend('bottomleft', leg=c('random','highest risk first','simultaneous instant'), lty = 1:3, lwd = 2, bty = 'n', title = 'order of cluster vaccination')
mtext('# simulations', 2, -1, outer = T)
mtext('proportion of district-level cases in trial population', 1, -1, outer = T)
graphics.off()
