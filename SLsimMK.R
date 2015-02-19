if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','SLSimsSWnoGLMM')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT')##,'RCT','FRCT','CRCT')
numEach <- 12*28

ves <- c(0, seq(.4, .9, by = .1))
pits <- c(.03, .05, .1)
parmsMat <- as.data.table(expand.grid(
    seed =  1:numEach
    , trial = tnms
    , ord = c('none','TU')
    , propInTrial = pits
    , sdLogIndiv = makeParms()$sdLogIndiv
    , delayUnit = 7#c(0,7)
    , vaccEff = 0#ves
    ))
parmsMat <- parmsMat[!(trial=='SWCT' & (delayUnit==0 | ord=='TU'))] ## SWCT must have delay and cannot be ordered
parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
nmtmp <- 'simSL-3-'
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 5
parmsMat$reordLag <- 14
parmsMat$nboot <- 200
nrow(parmsMat)

## fls <- list.files(batchdirnm, pattern=nmtmp)
## flsfull <- list.files(batchdirnm, pattern=nmtmp, full.names=T)
## sz <- unlist(sapply(fls, function(x) file.info(file.path(batchdirnm, x))['size']))
## fls <- fls[sz>6000]
## done <- gsub(nmtmp, '', fls)
## done <- as.numeric(gsub('.Rdata', '', done))
## length(done)

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

parmsMat[, length(nboot), propInTrial]
for(ii in 1:length(pits)) {
    parmsMatDo <- parmsMat[propInTrial==pits[ii]]
    ## parmsMatDo <- parmsMat[vaccEff==ves[ii]]
    sink(paste0('SLsims',pits[ii],'.txt'))
    for(ii in parmsMatDo$simNum) {
        cmd <- "R CMD BATCH '--no-restore --no-save --args"
        cmd <- addParm(cmd, parmsMatDo, ii)
        cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, formatC(ii, width=6, flag="0"),'.Rout')), sep='')
        cat(cmd)               # add command
        cat('\n')              # add new line
    }
    sink()
}

