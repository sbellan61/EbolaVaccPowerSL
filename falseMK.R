if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
 
batchdirnm <- file.path('BigResults','FalsePosFluct')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT','CRCT')
numEach <- 12*10

p2 <- simTrial(makeParms('RCT',small=F, ord='none', delayUnit = 0, clusSize=300, hazType = 'Phenom', weeklyDecay = .9, cvWeeklyDecay = .5, cvClus = 1.5, cvClusTime = 0.5, numClus = 20))

pits <- c(.025, .05, .075, .1)
parmsMat <- as.data.table(expand.grid(
    seed =  1:numEach
    , trial = tnms
    , ord = c('none','TU')
    , propInTrial = pits
    , sdLogIndiv = makeParms()$sdLogIndiv
    , delayUnit = 7
    , immunoDelay = 21
    , vaccEff = 0
    ))
parmsMat <- parmsMat[!(trial=='SWCT' & (delayUnit==0 | ord=='TU'))] ## SWCT must have delay and cannot be ordered
parmsMat <- parmsMat[!(delayUnit==0 & ord=='TU')] ## ordering is meaningless with simultaneous instant vacc
parmsMat <- parmsMat[ !(delayUnit==0 & trial=='FRCT')]  ## FRCT = RCT when delayUnit=0
parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
nmtmp <- 'simFP-'
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 17 ## 17*12 is ~ 2000 simulations each (2040 but we'll round)
parmsMat$reordLag <- 14
parmsMat$nboot <- 200
parmsMat$trialStartDate <- '2015-02-18'
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

parmsMat <- parmsMat[!(propInTrial!=.05 & immunoDelay!=21)]
parmsMat[, length(nboot), propInTrial]
parmsMat[, length(nboot), vaccEff]
parmsMat[, length(nboot), list(vaccEff,propInTrial)]
parmsMat[, length(nboot), list(immunoDelay,vaccEff,propInTrial)]
jbs <- NULL
immDs <- parmsMat[,unique(immunoDelay)]
jn <- 0
for(dd in 1:length(immDs)) {
    for(vv in 1:length(ves)) {
        for(pp in 1:length(pits)) {
            parmsMatDo <- parmsMat[vaccEff==ves[vv] & propInTrial==pits[pp] & immunoDelay==immDs[dd]]
            jn <- jn+1
            jbs <- rbind(jbs, data.table(jn=jn, vaccEff=ves[vv], propInTrial=pits[pp], immunoDelay=immDs[dd]))
            sink(paste0('SLsims',jn,'.txt'))
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
    }
}

jbs
