####################################################################################################
## Set up R CMD BATCH scripts for running on HPC cluster: Trial date sensitivity analysis
####################################################################################################
## Code base accompanying:
## 
## Bellan, SE, JRC Pulliam, CAB Pearson, DChampredon, SJ Fox, L Skrip, AP Galvani, M Gambhir, BA
## Lopman, TC Porco, LA Meyers, J Dushoff (2015). The statistical power and validity of Ebola
## vaccine trials in Sierra Leone: A simulation study of trial design and analysis. _Lancet
## Infectious Diseases_.
##
## Steve Bellan, March 2015
## License at bottom.
####################################################################################################

sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R'), source)
 
batchdirnm <- file.path('BigResults','initDateSens')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT')
numEach <- 12*10

trialSDseq <- paste0('2015-', c('01-15','02-01','02-15','03-01','03-15','04-01'))
ves <- .9
pits <- .05
parmsMat <- as.data.table(expand.grid( ## Matrix of parameter runs
    seed =  1:numEach
    , trial = tnms
    , propInTrial = pits
    , sdLogIndiv = makeParms()$sdLogIndiv
    , delayUnit = 7
    , immunoDelay = 21 
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
nrow(parmsMat)

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
            sink(paste0('initDateSens',jn,'.txt')) ## file to be submitted to the cluster
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

####################################################################################################
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/sbellan61/EbolaVaccPowerSL
####################################################################################################
