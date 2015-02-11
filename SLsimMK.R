if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSims3')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT','CRCT')
numEach <- 4

parmsMat <- as.data.table(expand.grid(
    seed =  1:numEach
    , trial = tnms
    , propInTrial = c(.03, .05, .1)
    , sdLogIndiv = makeParms()$sdLogIndiv
    , delayUnit = c(0,7)
    , vaccEff = c(0, seq(.4, .95, by = .05))
    ))
parmsMat <- parmsMat[!(trial=='SWCT' & delayUnit==0)] ## SWCT must have delay
parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
nmtmp <- 'simSL-3-'
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 300
parmsMat$ord <- "none"
parmsMat[trial %in% c('RCT','FRCT','CRCT'), ord := 'TU']
parmsMat$reordLag <- 14
nrow(parmsMat)

## fls <- list.files(batchdirnm, pattern=nmtmp)
## sz <- unlist(sapply(fls, function(x) file.info(file.path(batchdirnm, x))['size']))
## fls <- fls[sz>0]
## done <- gsub(nmtmp, '', fls)
## done <- as.numeric(gsub('.Rdata', '', done))
## length(done)

## parmsMatDo <- parmsMat[!simNum %in% done]
## parmsMatDo <- parmsMat[trial %in% c('RCT','FRCT')]
## nrow(parmsMatDo)

addParm <- function(x, parmsMat,ii) {
    for(pp in 1:length(parmsMat)) {
        tempP <- as.data.frame(parmsMat)[,pp]
        isch <- !is.numeric(tempP[1])
        parmAdd <- tempP[parmsMat$simNum==ii]
        addStrg <- paste0(" ", names(parmsMat)[pp], "=", "\""[isch], parmAdd, "\""[isch])
        x <- paste0(x, addStrg)
    }
    return(x)
}

sink('SLsims.txt')
for(ii in parmsMatDo$simNum) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMatDo, ii)
    cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, formatC(ii, width=6, flag="0"),'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
