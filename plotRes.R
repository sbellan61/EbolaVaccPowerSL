if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R'), source)
## Analyze stopping time results
batchdirnm <- file.path('Results','StpSims')
fls <- list.files(batchdirnm, pattern='Rdata', full.names=T)
load(fls[1])
stps

pdf(file.path('Figures','stopping time results.pdf'))
par(mfrow=c(4,1), mar = c(4,4,1,0))
xlim <- c(0,200)
for(typ in c('swct','rct','frct','crct')) {
    tempout <- as.data.table(get(paste0(typ,'Sims')))
    hist(tempout[,stopDay], breaks = seq(0, 10^4, by = 10), col = 'black', xlab ='stopping time', main = typ, xlim = xlim)
    abline(v=mean(tempout[,stopDay]), col = 'red')
}
graphics.off()
