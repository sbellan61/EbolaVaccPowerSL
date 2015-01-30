batchdirnm <- file.path('Results','StpSims')
routdirnm <- file.path('Results','StpSims','Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT','CRCT')
trials <- rep(tnms, each = 12)
nsims <- 100

sink('stpSims.txt')
todo <- 1:length(trials)
for(ii in todo) {
    cmd <- paste("R CMD BATCH '--no-restore --no-save --args batchdirnm=\"", batchdirnm, "\" simNum=", ii, " nsims=", nsims, " trial=\"", trials[ii],
                 "\" ' stopTimes.R ", file.path(batchdirnm,'Routs', paste0('stpSim', ii,'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
