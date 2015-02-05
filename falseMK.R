if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
batchdirnm <- file.path('Results','FalsePosSims')
routdirnm <- file.path('Results','FalsePosSims','Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT','CRCT')
numEach <- 12

parmsMat <- data.frame(
    nsims = 200
  , seed =  rep(1:numEach, length(tnms))
  , trial = rep(tnms, each = numEach)
  , batchdirnm = batchdirnm
  , varClus = 0
  , weeklyDecay = 1
  , weeklyDecayVar = 0
  , vaccEff = 0
)
parmsMat$simNum <- 1:nrow(parmsMat)

addParm <- function(x, parmsMat,ii) {
    for(pp in 1:length(parmsMat)) {
        tempP <- parmsMat[,pp]
        isch <- !is.numeric(tempP[1])
        parmAdd <- tempP[ii]
        addStrg <- paste0(" ", names(parmsMat)[pp], "=", "\""[isch], parmAdd, "\""[isch])
        x <- paste0(x, addStrg)
    }
    return(x)
}

sink('falsePosSims.txt')
for(ii in parmsMat$simNum) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMat, ii)
    cmd <- paste0(cmd, " ' falsePos.R ", file.path(batchdirnm,'Routs', paste0('falsePosSim', ii,'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
