batchdirnm <- file.path('Results','FalsePosSims')
routdirnm <- file.path('Results','FalsePosSims','Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT','CRCT')
numEach <- 12
trials <- rep(tnms, each = numEach)
seeds <- rep(1:numEach, length(tnms))
nsims <- 100

sink('falsePosSims.txt')
todo <- 1:length(trials)
for(ii in todo) {
    cmd <- paste("R CMD BATCH '--no-restore --no-save --args batchdirnm=\"", batchdirnm, "\" simNum=", ii, " nsims=", nsims, " trial=\"", trials[ii],
                 "\" ' stopTimes.R ", file.path(batchdirnm,'Routs', paste0('falsePosSim', ii,'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
