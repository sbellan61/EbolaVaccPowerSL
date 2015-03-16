if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
 
batchdirnm <- file.path('BigResults','initDateSens')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT')#,'CRCT')
numEach <- 12*10

trialSDseq <- paste0('2015-', c('01-15','02-01','02-15','03-01','03-15','04-01'))
ves <- .9
pits <- .05
parmsMat <- as.data.table(expand.grid(
    seed =  1:numEach
    , trial = tnms
    ##    , ord = c('none','TU')
    , propInTrial = pits
    , sdLogIndiv = makeParms()$sdLogIndiv
    , delayUnit = 7
    , immunoDelay = 21 ## c(5, 21)
    , vaccEff = ves
    , trialStartDate = trialSDseq
    ))
parmsMat$remProtDel <- T
parmsMat$remStartFin <- F
parmsMat$ord <- 'TU'
parmsMat[trial=='SWCT', ord := 'none']
parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
nmtmp <- 'simSL-initD-remPD-'
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 17 ## 17*12 is ~ 2000 simulations each (2040 but we'll round)
parmsMat$reordLag <- 14
parmsMat$nboot <- 200
## parmsMat$trialStartDate <- '2015-02-18'
nrow(parmsMat)

## fls <- list.files(batchdirnm, pattern=nmtmp)
## flsfull <- list.files(batchdirnm, pattern=nmtmp, full.names=T)
## sz <- unlist(sapply(fls, function(x) file.info(file.path(batchdirnm, x))['size']))
## fls <- fls[sz>6000]
## done <- gsub(nmtmp, '', fls)
## done <- as.numeric(gsub('.Rdata', '', done))
## parmsMat <- parmsMat[simNum[!simNum %in% done]]
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

parmsMat <- parmsMat[trial=='SWCT']
parmsMat[, length(nboot), trial]
jbs <- NULL
jn <- 0
for(tt in 1:length(tnms)) {
            parmsMatDo <- parmsMat[trial==tnms[tt]]
            jn <- jn+1
            jbs <- rbind(jbs, data.table(jn=jn, trial=tnms[tt]))
            sink(paste0('initDateSens',jn,'.txt'))
            for(ii in parmsMatDo$simNum) {
                cmd <- "R CMD BATCH '--no-restore --no-save --args"
                cmd <- addParm(cmd, parmsMatDo, ii)
                cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, sprintf("%06d", ii),'.Rout')), 
                              sep='')
                cat(cmd)               # add command
                cat('\n')              # add new line
            }
            sink()
}
jbs
