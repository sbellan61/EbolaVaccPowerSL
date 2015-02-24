if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid)
##load(file=file.path('BigResults','powFin.Rdata'))
percent <- function(x) paste0(formatC(x*100), '%')
labs <- c('','log')

thing <- 'initDateSens'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))

pf[delayUnit==0, ord:='simultaneous instant']
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.025 & mod=='coxME']
setnames(pf,'ord','order')
pf$design <- pf$trial
levels(pf$design)[levels(pf$design) == 'SWCT'] <- 'SWT'
levels(pf$order)[2] <- 'time-updated'
pf[, immunoDelay:=as.numeric(levels(immunoDelay)[immunoDelay])]
pf[, pit:=factor(paste0(propInTrial*100,'%'))]
pf[, pit:=factor(pit, levels = c('2.5%','5%','7.5%','10%'), ordered = T)]
## levels(pf$pit) <- c('2.5%','5%','7.5%','10%')

baseMods <- c('Cox PH Frailty'
              , 'Poisson GLM\n no cluster effects'
              , 'Poisson GLM \nwith fixed effects by cluster')

pf$model <- pf$mod
levels(pf$model) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))


####################################################################################################
## Figure 4 - Power
subs <- pf[, immunoDelay==21 & ((design == 'SWT' & mod=='relabCoxME') | (design == 'RCT' & order=='time-updated' & mod =='coxME'))]
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[subs], aes(trialStartDate, vaccGoodNAR, colour=design, linetype=order)) + thsb +
            scale_x_date(labels = date_format("%b-%d"), breaks = pf[,unique(trialStartDate)], minor_breaks=NULL) +  
                          xlab('trial start date') + ylab('power to detect effective vaccine') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
                                  geom_line(size=1) #+ facet_wrap(~pit, scales = "free_y",nrow=1) + 
#                                      ggtitle('expected % of district-level cases in trial population')
p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
ggsave(paste0('Figures/Fig 6 - Power by start date SL.png'), p.tmp, w = 6, h = 4)

