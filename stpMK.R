batchdirnm <- file.path('Results','Sims')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(batchdirnm)) dir.create(file.path(batchdirnm,'Routs'))
tnms <- c('SWCT','RCT','FRCT','CRCT')
trials <- rep(tnms, each = 3)
nsims <- 10

sink('stpSims.txt')
todo <- 1:length(trials)
for(ii in todo) {
    cmd <- paste("R CMD BATCH '--no-restore --no-save --args simnum=", ii, " nsims=", nsims, " trial=", trials[ii],
                 " stopTimes.R ", file.path(batchdirnm,'Routs', paste0('stpSim', ii,'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
