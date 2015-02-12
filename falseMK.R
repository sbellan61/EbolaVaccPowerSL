if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('BigResults','FalsePosSims2')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','FRCT') #,'FRCT','CRCT')
numEach <- 8

parmsMat <- as.data.table(expand.grid(
    seed =  5:(5+numEach)
    , trial = tnms
    , sdLogIndiv = makeParms()$sdLogIndiv
    , cvClus = c(0, 1, 1.5)
    , weeklyDecay = .9
    , cvWeeklyDecay = c(0, .5, 1, 1.5)
    , vaccEff = c(0, seq(.4, .95, by = .05))
    ))
parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
parmsMat$hazSL <- F
nmtmp <- 'simFP-big2'
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 600
nrow(parmsMat)

## fls <- list.files(batchdirnm, pattern=nmtmp)
## done <- gsub(nmtmp, '', fls)
## done <- as.numeric(gsub('.Rdata', '', done))
## length(done)

parmsMatDo <- parmsMat ##[!simNum %in% done]
parmsMatDo <- parmsMat[vaccEff==0]
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
    cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, formatC(ii, width=6, flag="0"),'.Rout')), sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
