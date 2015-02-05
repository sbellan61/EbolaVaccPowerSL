if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('Results','FalsePosSims')
routdirnm <- file.path('Results','FalsePosSims','Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT','CRCT')
numEach <- 12

parmsMat <- as.data.table(expand.grid(
    seed =  1:numEach
    , trial = tnms
    , sdLogIndiv = c(0,makeParms()$sdLogIndiv)
    , varClus = c(0, makeParms()$varClus) 
    , weeklyDecay = c(1, makeParms()$weeklyDecay)
    , weeklyDecayVar = c(0, makeParms()$weeklyDecayVar)
    , vaccEff = c(0, .6, .7, .8)
    ))
parmsMat <- parmsMat[! ( weeklyDecay==1 & weeklyDecayVar!=0 )] ## not interested in variance around stability
parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
parmsMat$saveNm <- 'simFP2-'
parmsMat$nsims <- 100

parmsMatDo <- parmsMat[sdLogIndiv==1 & varClus!=0 & weeklyDecay!=1 & weeklyDecayVar!=0]
nrow(parmsMatDo)

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

sink('falsePosSims.txt')
for(ii in parmsMatDo$simNum) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMatDo, ii)
    cmd <- paste0(cmd, " ' falsePos.R ", file.path(batchdirnm,'Routs', paste0('falsePosSim', ii,'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
