if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)
##load(file=file.path('BigResults','powFin.Rdata'))
percent <- function(x) paste0(formatC(x*100), '%')
labs <- c('','log')

load(file=file.path('Results','powFin_All.Rdata'))

pf[vaccEff==.5 & trial=='RCT' & propInTrial==.025 & mod=='CoxME']
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
group.colors <- c(RCT = "#333BFF", FRCT = "#CC6600", SWT ="#9633FF")
group.colors[c(1,3,2)] <- gg_color_hue(3)
group.colors['SWCT'] <- 'orange'
pf$trial <- factor(pf$trial, levels=levels(pf$trial)[c(2,1,3)])
pf[, biasNAR:=biasNAR/vaccEff]
levels(pf$order)[1] <- 'random'
levels(pf$order)[2] <- 'risk-prioritized'

## Set up line types for ease
pf$analysis <- 1
pf[trial=='SWCT' & swctPT=='remPD', analysis:=1]
pf[trial=='SWCT' & swctPT=='all', analysis:=2]
pf[trial=='SWCT' & swctPT=='remPD_SF', analysis:=3]
pf$analysis <- factor(pf$analysis, levels = 1:3, labels = c(
                                                     'SWCT including \nearly/late person-time \n& excluding protective delay',
                                                     'SWCT including protective delay \n& early/late person-time',
                                                     'SWCT excluding \nboth protective delay \n& early/late person-time'))
## Relabel models
pf$modLab2 <- pf$mod
levels(pf$modLab2)[levels(pf$modLab2) %in% c('CoxME','bootCoxME','relabCoxME')] <- c('(A) CoxPH', '(B) Bootstrap', '(C) Permutation')
## Only include SWCT remPD_SF pt
pf$mainAn <- pf[, trial!='SWCT' | (trial=='SWCT' & swctPT=='remPD')]

####################################################################################################
## them for ms
thax <- element_text(colour = 'black', size = 8)
thsb <- theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
              axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
              axis.line = element_line(), axis.ticks = element_line(color='black'),
              panel.margin = unit(1, "lines"), legend.key.height=unit(2,"line")
              , strip.background = element_rect(fill = NA)
              ,legend.position = 'right'
              , axis.line = element_blank()
              ,panel.grid.major = element_blank()
              , panel.grid.minor = element_blank()
              ,panel.border = element_blank()
              ,panel.background = element_blank()
              , legend.background =  element_blank()
              , legend.key =  element_blank()
              , legend.key.width=unit(2,"line")
              ## ,legend.justification=c(1,0), legend.position=c(1,0)
              )
theme_set(theme_grey(base_size = 12))
##thsb <- thsb + theme_bw()#

####################################################################################################
## Figure 4 - Type I errors
subs <- pf[, mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=trial, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.025,.1), minor_breaks=NULL, breaks = c(.025,.05,.075,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('False Positive Rate') + 
    #scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~modLab2, scales = "free_y") + scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
p.tmpunt
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Fig 4 - Type I by analysis',typ), p.tmpunt, w = 8.5, h = 3.5)

subs <- pf[, mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','GLMFclus', 'bootCoxME','bootGLMFclus', 'relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
pf[vaccEff==0 & immunoDelay==21 & subs & pit=='10%', list(trial, stoppedNAR,mod)]

####################################################################################################
## Comparison of SWCT pt stuff
subs <- pf[, trial %in% c('SWCT') & vaccEff==0 & mod %in% c('CoxME','GLMFclus','bootCoxME','relabCoxME') & immunoDelay==21]
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = percent, limits=c(0,.3), minor_breaks=NULL, breaks = c(0,.1,.2,.3)) +  
    xlab('% of district-level cases in trial population') + ylab('False Positive Rate') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) + theme(legend.key.height=unit(3,"line")) + 
    geom_line(size=1, colour = group.colors['SWCT']) + facet_wrap(~modLab2, scales = "free_y", nrow = 1) + scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.17))
p.tmpunt
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/type I by SWCT person-time',typ), p.tmpunt, w = 12, h = 3.5)

####################################################################################################
## Figure 4B - Type I errors panels by trial type
thax <- element_text(colour = 'black', size = 8)
subs <- pf[, mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & mod %in% c('CoxME','bootCoxME','relabCoxME') & immunoDelay==21]
subs <- subs & pf[,!(trial=='RCT' & grepl('boot',mod))]
temp.cols <- c('Cox PH Frailty'='orange', 'permutation test over\nCox PH Frailty'='purple',
               'bootstrap over\nCox PH Frailty' = 'dodger blue')
p.tmp <- ggplot(pf[subs], 
                aes(propInTrial, stoppedNAR, colour=model, linetype=order)) + thsb +
    scale_x_continuous(labels = percent, limits=c(.025,.1), minor_breaks=NULL, breaks = c(.025,.05,.075,.1)) +  
    xlab('% of district-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_hline(yintercept=.05, color='dark gray', size = 1) +
    geom_line(size=1) + facet_wrap(~trial, scales = "free_y") + scale_color_manual(values=temp.cols)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,.15))
p.tmpunt
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Type I by trial SL', typ), p.tmpunt, w = 6.5, h = 3.5)

####################################################################################################
## Figure 5 - Power
subs <- pf[,  mainAn==T & vaccEff>0 &  immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT','FRCT') & mod =='CoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=trial, linetype=order)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit, scales = "free_y",nrow=1) + 
    ggtitle('expected % of district-level cases in trial population') + 
    scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/Fig 5 - Power SL', typ), p.tmpunt, w = 9, h = 4)

tmpS <- arrange(pf[subs & trial=='SWCT', list(trial, vaccEff, vaccGoodNAR, order,pit)],pit)
tmpR <- arrange(pf[subs & trial=='RCT' & order=='risk-prioritized', list(trial, vaccEff, vaccGoodNAR, order,pit)], pit)
tmpR$vaccGoodNAR/tmpS$vaccGoodNAR

####################################################################################################
##  Power by immunedelay
for(ii in 1:2) {
    if(ii==1) pdf('Figures/Fig 5 - Power by seroconversion delay SL.pdf', w = 4, h=5)
    if(ii==2) png('Figures/Fig 5 - Power by seroconversion delay SL.png', w = 4, h=5, res = 300, units='in')
##    layout(matrix(1:2,ncol=2), w = c(1,.5))
    par(mar = c(5,4,0,0))
    subs <- pf[, propInTrial==.05 & vaccEff == .9 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial=='RCT' & order=='risk-prioritized' & mod =='CoxME'))]
    pf[subs]
    plot(0,0, type = 'n', xlim = c(5, 21), ylim = c(0,1), bty='n', las = 1, axes=F,
         xlab ='delay until protective efficacy', ylab='power to detect effective vaccine')
    axis(1, c(5,21))
    axis(2, seq(0,1,l=5),las=1)
    with(pf[subs & trial=='RCT'], lines(immunoDelay, vaccGoodNAR, col = group.colors['RCT'], lty = 2, lwd = 3))
    with(pf[subs & trial=='SWCT' & swctPT=='remPD'], lines(immunoDelay, vaccGoodNAR, col = group.colors['SWCT'], lty = 1, lwd = 3))
    with(pf[subs & trial=='SWCT' & swctPT=='all'], lines(immunoDelay, vaccGoodNAR, col = group.colors['SWCT'], lty = 2, lwd = 3))
    with(pf[subs & trial=='SWCT' & swctPT=='remPD_SF'], lines(immunoDelay, vaccGoodNAR, col = group.colors['SWCT'], lty = 3, lwd = 3))
    ## par(mar=rep(0,4))
    ## plot.new()
    ## legend('top', leg = c('risk-prioritized RCT', levels(pf$analysis)),
    ##        col = c(group.colors['RCT'], rep(group.colors['SWCT'],3)), lty = c(2,1:3), lwd = 3, bty = 'n', ncol = 2, cex = .8)
    graphics.off()
}

####################################################################################################
##  Power by SWCT design
subs <- pf[, vaccEff>0 & propInTrial <= .1 & immunoDelay==21 & (trial == 'SWCT' & mod %in% c('CoxME','GLMFclus','bootCoxME','relabCoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(.5,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_line(size=1, colour=group.colors['SWCT']) + facet_wrap(mod~pit, scales = "free",nrow=4) +
 theme(legend.key.height=unit(3,"line")) + 
    ggtitle('expected % of district-level cases in trial population') + 
    scale_color_manual(values=group.colors)
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/Power by SWCT pt', typ), p.tmpunt, w = 9, h = 8)

####################################################################################################
##  Prob rej null by SWCT design
subs <- pf[, immunoDelay==21 & propInTrial <= .1 & (trial == 'SWCT' & mod %in% c('CoxME','bootCoxME','relabCoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, linetype=analysis)) + thsb +
    scale_x_continuous(labels = formatC, limits=c(0,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('power to detect effective vaccine') + 
    scale_linetype_manual(breaks=levels(pf$analysis), values=1:3) +
    geom_line(size=1, colour=group.colors['SWCT']) + facet_wrap(mod~pit, scales = "free",nrow=3) +
    theme(legend.key.height=unit(3,"line")) + 
    ggtitle('expected % of district-level cases in trial population') + 
    geom_hline(yintercept=.025, color='dark gray', size = 1) + 
    scale_color_manual(values=group.colors)
p.tmplog <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25, .5, 1), limits = c(.005,1), minor_breaks=NULL) 
p.tmplog
for(typ in c('.png','.pdf'))    ggsave(paste0('Figures/CoxPH Prob rej null by SWCT pt log', typ), p.tmplog, w = 9, h = 8)

####################################################################################################
## Figure SX - Power by model
pf$modelLab <- pf$model
pf$class <- 'parametric'
pf[grepl('boot',mod), class:='bootstrap']
pf[grepl('relab',mod), class:='permutation']
pf$class <- factor(pf$class, levels = c("parametric", "bootstrap", "permutation"))
pf$model <- gsub('boot','',pf$mod)
pf$model <- factor(gsub('relab','',pf$model))
levels(pf$model)[1:2] <- c('CoxPH','Poisson regression')
subs <- pf[, mainAn==T & immunoDelay==21 & (trial == 'SWCT' | (trial=='RCT' & order=='risk-prioritized' & class!='bootstrap'))]
subs <- subs & pf[,model %in% c('CoxPH','Poisson regression')]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=class, linetype=model)) +
    thsb +
    scale_x_continuous(labels = formatC, limits=c(0,1), breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('probability of rejecting null that vaccine is not effective') + 
    scale_linetype_manual(values=1:3) + ## scale_color_manual(values=group.colors) +
    geom_hline(yintercept=.025, color='dark gray', size = 1) + 
    scale_color_manual(values=c(parametric='orange',bootstrap='dodger blue',permutation='purple')) + 
    geom_line(size=.9)  + facet_wrap(trial ~ pit, scales = "free_y", nrow = 2) #,nrow=1) + 
#p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
p.tmpunt <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1))
p.tmplog <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .25, .5, 1), limits = c(.005,1), minor_breaks=NULL) 
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Power by model class',typ), p.tmpunt, w = 8, h = 6)
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Power by model class log',typ), p.tmplog, w = 8, h = 6)

##################################################
## Cases
subs <- pf[, mainAn==T & immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT','FRCT') & mod =='CoxME'))]
thax <- element_text(colour = 'black', size = 8)
p.tmp <- ggplot(pf[subs], aes(vaccEff, caseTot, colour=trial, linetype=order)) + thsb + 
    scale_x_continuous(labels = formatC, limits=c(0,.9), breaks = pf[,unique(vaccEff)],minor_breaks=NULL) +  
    xlab('vaccine efficacy') + ylab('# EVD cases') +
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~pit, scale='free_y', nrow=1) + 
    ggtitle('% of district-level cases in trial population')  + scale_color_manual(values=group.colors) + 
    scale_y_continuous(labels = formatC, limits=c(0,100)) 
## subs <- pf[, swctPT=='remPD_SF' & immunoDelay==21 & (trial == 'SWCT' & mod=='relabCoxME')]
## p.tmp <- p.tmp + geom_line(aes(vaccEff, caseTot), data = pf[subs], linetype = 4, col = 'orange', size = 1.3)
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Cases SL ALL',typ), p.tmp, w = 9, h = 3.5)

subs <- pf[,vaccEff==0 & (delayUnit==0 | (trial=='SWCT')) & mod =='CoxME']
pf[subs, list(caseTot,caseC,caseV,pit,trial,delayUnit,vaccEff)]

####################################################################################################
## abstract #'s and tables
pf <- arrange(pf, immunoDelay, trial, propInTrial, mod)

subs <- pf[,mainAn==T & immunoDelay==21 & (trial=='SWCT' & mod=='relabCoxME') | (trial=='RCT' & ord=='TU' & mod=='CoxME')]
pf[subs & vaccEff==.9, list(trial, ord, delayUnit, mod, vaccGoodNAR, propInTrial,vaccEff, immunoDelay)]

pf[mainAn==T  & immunoDelay==21 & mod %in% c('CoxME') & vaccEff==0 &  (trial=='SWCT' | ( trial=='RCT' & ord=='TU')),
list(trial, ord, delayUnit, mod, stoppedNAR, propInTrial,vaccEff)]

pf[mainAn==T  & mod %in% c('CoxME') & vaccEff>.5 &  (trial=='SWCT' | ( trial=='FRCT' & ord=='TU')),
list(trial, ord, mod, vaccGoodNAR,cvr, propInTrial,vaccEff)]

## Chosen models
subs <- pf[, mainAn==T  & vaccEff %in% c(0,.9) & immunoDelay==21 & #propInTrial == .05 & 
           ((mod %in% c('CoxME','relabCoxME') & trial=='SWCT' ) |
           (mod %in% c('CoxME') & trial=='RCT' & ord=='TU' ))]
tab1 <- pf[subs, list(trial, mod, pit, cvr, biasNAR, vaccEff)]
tab1
tab1[, list(trial, cvr, biasNAR)]
tab1[trial=='SWCT', cvr:= cvr[!is.na(cvr)], list(pit,vaccEff)]
tab1 <- tab1[!(trial=='SWCT' & mod=='CoxME')]
tab1[, c('cvr','biasNAR') := list(signif(cvr,2), signif(biasNAR,2))]
tab1

tab1 <- data.table(tab1[trial=='RCT',list(pit, cvr,biasNAR)], tab1[trial=='SWCT',list(cvr,biasNAR)])
tab1
write.csv(tab1, file='Results/biasCoverage.csv')

#################################################
## Power
for(jj in 1:2) {
    pdf(paste0('Figures/',labs[jj], 'power SL ALL.pdf'), w = 8, h = 4) ##, units = 'in', res = 200)
    for(modTmp in levels(pf$model)) {
        subs <- pf[, mainAn==T & model==modTmp & immunoDelay==21 & !(trial%in%c('FRCT','RCT') & grepl('boot',mod))]
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[subs], aes(vaccEff, vaccGoodNAR, colour=trial, linetype=order)) + 
            scale_x_continuous(labels = formatC, limits=c(0,.9),  breaks = pf[,unique(vaccEff)], minor_breaks=NULL) +  
                theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
                      axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
                      panel.margin = unit(2, "lines"))  + scale_color_manual(values=group.colors) + 
                          xlab('vaccine efficacy') + ylab('probability of rejecting null') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                  geom_line(size=1) + facet_wrap(~pit, nrow=1) + 
                                      ggtitle('power by expected % of district-level cases in trial population') +
                                          geom_hline(yintercept=.025, color='black', linetype='dotted')
        subs <- pf[, swctPT=='remPD_SF' & model==modTmp & immunoDelay==21 & trial=='SWCT']
        p.tmp <- p.tmp + geom_line(aes(vaccEff, vaccGoodNAR), colour=group.colors['SWCT'], linetype=4, data =pf[subs], size = 1)
        if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), minor_breaks = seq(0,1,.05))
        if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, limits=c(0.0005,.96), breaks = c(.01,.025,.05,.1,.2,.5,.8))
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
subs <- pf[,mainAn==T & trial %in% c('SWCT','RCT') & vaccEff==0 & immunoDelay==21 & !(trial=='RCT' & grepl('boot',mod))]
p.tmp <- ggplot(pf[subs], aes(propInTrial, stoppedNAR, colour=trial, linetype=order)) + 
    scale_x_continuous(labels = percent, limits=c(.025,.1), breaks = c(.025,.05,.075,.1)) +  
##    scale_y_continuous(labels = formatC, limits=c(0,.3)) + #, minor_breaks = seq(0,1,.05)) +
    theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
          axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
          panel.margin = unit(2, "lines")) + scale_color_manual(values=group.colors) +
    xlab('% of distrinct-level cases in trial population') + ylab('Type I Error Rate') + 
    scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
    geom_line(size=1) + facet_wrap(~model) + 
geom_hline(yintercept=.05, color='black', linetype='dotted', size = 1.2)
       subs <- pf[,swctPT=='remPD_SF' & trial =='SWCT' & vaccEff==0 & immunoDelay==21 & !(trial=='RCT' & grepl('boot',mod))]
        p.tmp <- p.tmp + geom_line(aes(propInTrial, stoppedNAR), colour=group.colors['SWCT'], linetype=4, data =pf[subs], size = 1)
if(jj==1) p.tmp <- p.tmp + scale_y_continuous(labels = formatC)#, limits=c(0,.3))
if(jj==2) p.tmp <- p.tmp + scale_y_log10(labels = formatC, breaks = c(.01, .025, .05, .1, .2, .3, .4))#,  limits=c(.01,.3))
print(p.tmp)
graphics.off()
}

####################################################################################################
## Show power by # of cases
subs <- pf[,  mainAn==T & vaccEff==.9 & immunoDelay==21 & ((trial == 'SWCT' & mod=='relabCoxME') | (trial %in% c('RCT') & mod =='CoxME'))]
pf$type <- pf[, paste(order,trial)]
pf[grepl('SWCT',type), type:='SWCT']
group.colors2 <- c('purple','blue', "#F8766D", "orange")
names(group.colors2) <- pf[subs,unique(type)]
p.tmp <- ggplot(pf[subs], aes(caseTot, vaccGoodNAR, colour=type)) + thsb +
    xlab('# of cases in trial population') + ylab('power to detect effective vaccine') + 
    geom_line(size=1.5)  +
    geom_point(aes(caseTot, vaccGoodNAR, colour=type, shape=pit), size = 5) +
    scale_color_manual(values=group.colors2) + scale_shape_manual(values=c(15:18)) + 
    scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) 
subs <- pf[,  swctPT=='remPD_SF' & vaccEff==.9 & immunoDelay==21 & (trial == 'SWCT' & mod=='relabCoxME')]
p.tmp <- p.tmp + geom_point(aes(caseTot, vaccGoodNAR, shape=pit), data = pf[subs], size = 5) +
    geom_line(aes(caseTot, vaccGoodNAR), data = pf[subs], linetype=4, size = 1.3)
for(typ in c('.png','.pdf')) ggsave(paste0('Figures/Fig SX - Power by cases',typ), p.tmp, w = 8, h = 4)


