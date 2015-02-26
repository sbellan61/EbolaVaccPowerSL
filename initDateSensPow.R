if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(data.table); library(ggplot2); library(dplyr); library(grid); library(scales)
##load(file=file.path('BigResults','powFin.Rdata'))
percent <- function(x) paste0(formatC(x*100), '%')
labs <- c('','log')

thing <- 'initDateSens'
load(file=file.path('Results',paste0('powFin_',thing,'.Rdata')))
pf[vaccEff==.5 & trial=='RCT' & propInTrial==.025 & mod=='CoxME']
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
group.colors <- c(RCT = "#333BFF", FRCT = "#CC6600", SWT ="#9633FF")
group.colors[c(1,3,2)] <- gg_color_hue(3)
group.colors['SWT'] <- 'orange'
pf$design <- factor(pf$design, levels=levels(pf$design)[c(2,1,3)])
pf[, biasNAR:=biasNAR/vaccEff]

####################################################################################################
## them for ms
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
              ## ,legend.justification=c(1,0), legend.position=c(1,0)
              )
theme_set(theme_grey(base_size = 12))
##thsb <- thsb + theme_bw()#

####################################################################################################
## Figure 4 - Power
subs <- pf[, immunoDelay==21 & ((design == 'SWT' & mod=='relabCoxME') | (design == 'RCT' & order=='time-updated' & mod =='CoxME'))]
        thax <- element_text(colour = 'black', size = 8)
        p.tmp <- ggplot(pf[subs], aes(trialStartDate, vaccGoodNAR, colour=design, linetype=order)) + thsb +
            scale_x_date(labels = date_format("%b-%d"), breaks = pf[,unique(trialStartDate)], minor_breaks=NULL) +  
                          xlab('trial start date') + ylab('power') + 
                              scale_linetype_manual(breaks=levels(pf$order), values=1:3) +
 theme(axis.text.x = element_text(angle=90)) + scale_color_manual(values=group.colors) +
     geom_segment(aes(x = as.Date('2015-02-18'), y =1, xend = as.Date('2015-02-18'), yend = .77), 
                  color='black',arrow = arrow(length = unit(0.2, "cm"))) +
     geom_segment(aes(x = as.Date('2015-03-18'), y =.7, xend = as.Date('2015-03-18'), yend = .5), 
                  color='black',arrow = arrow(length = unit(0.2, "cm"))) +
                                  geom_line(size=1) #+ facet_wrap(~pit, scales = "free_y",nrow=1) + 
#                                      ggtitle('expected % of district-level cases in trial population')
p.tmp <- p.tmp + scale_y_continuous(labels = formatC, limits=c(0,1), breaks=seq(0,1,by=.1), minor_breaks = NULL) #seq(0,1,.05))
ggsave(paste0('Figures/Fig 6 - Power by start date SL.png'), p.tmp, w = 5, h = 3)


subs <- pf[, vaccEff>.8 & immunoDelay==21 & propInTrial == .05 & 
           ((mod %in% c('relabCoxME') &  (trial=='SWCT' & ord=='none') ) |
           (mod %in% c('CoxME') & trial=='RCT' & ord=='TU' ))]
tab1 <- pf[subs, list(vaccGoodNAR), list(trial, pit,trialStartDate)]
tab1 <- dcast.data.table(tab1, trial ~ trialStartDate)
tab2 <- t(tab1[,-1,with=F])
tab2
colnames(tab2) <- tab1[,trial]
tab2[-1,]/tab2[-nrow(tab2),]

1-tab2[5,]/tab2[3,]

