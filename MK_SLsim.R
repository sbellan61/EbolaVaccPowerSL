####################################################################################################
## Set up R CMD BATCH scripts for running on HPC cluster: Main analysis.
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
 
batchdirnm <- file.path('BigResults','All')
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT')
numEach <- 12*10

ves <- c(0, seq(.5, .9, by = .2))
pits <- c(.025, .05, .075, .1, .15,.2,.3)
parmsMat <- as.data.table(expand.grid(
    seed =  1:numEach
    , trial = tnms
    , ord = c('none','TU')
    , propInTrial = pits
    , sdLogIndiv = makeParms()$sdLogIndiv
    , delayUnit = c(0,7)
    , immunoDelay = c(5, 21)
    , vaccEff = ves
    , remStartFin = c(F,T)
    , remProtDel = c(F,T)    
    ))
parmsMat <- parmsMat[!(trial=='SWCT' & (delayUnit==0 | ord=='TU'))] ## SWCT must have delay and cannot be ordered
parmsMat <- parmsMat[!(delayUnit==0 & ord=='TU')] ## ordering is meaningless with simultaneous instant vacc
parmsMat <- parmsMat[ !(delayUnit==0 & trial=='FRCT')]  ## FRCT = RCT when delayUnit=0
parmsMat <- parmsMat[ !(propInTrial>.1 & trial!='SWCT')] ## only for SWCT sensitivity analysis
parmsMat <- parmsMat[!(remProtDel==F & remStartFin==T)] ## not an option considered for SWCT persontime
parmsMat <- parmsMat[!(trial!='SWCT' & (remProtDel!=F | remStartFin!=F))] ## parameters only apply to SWCT
parmsMat[trial!='SWCT', remProtDel := NA]
parmsMat[trial!='SWCT', remStartFin := NA]
parmsMat <- parmsMat[!(immunoDelay==5 & propInTrial!=0.05)] ## only doing this sensitivity analysis for propInTrial of 5%

parmsMat$simNum <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
nmtmp <- 'simSL-main-'
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 17 ## 17*12 is ~ 2000 simulations each (2040 but we'll round)
parmsMat$reordLag <- 14
parmsMat$nboot <- 200
parmsMat$trialStartDate <- '2015-02-18'
nrow(parmsMat)

nrow(parmsMat)*17 ## 399,840 simulations

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
parmsMat <- parmsMat[!(propInTrial!=.05 & immunoDelay!=21)]
parmsMat[, length(nboot), propInTrial]
parmsMat[, length(nboot), vaccEff]
parmsMat[, length(nboot), list(vaccEff,propInTrial)]
parmsMat[, length(nboot), list(immunoDelay,vaccEff,propInTrial)]
nrow(parmsMat)
jbs <- NULL
immDs <- parmsMat[,unique(immunoDelay)]
jn <- 0

## for(dd in 1:length(immDs)) {
##for(vv in 1:length(ves)) {
        ## for(pp in 1:length(pits)) {
            parmsMatDo <- parmsMat#[vaccEff==ves[vv]]# & propInTrial==pits[pp] & immunoDelay==immDs[dd]]            
            jn <- jn+1
##            jbs <- rbind(jbs, data.table(jn=jn, vaccEff=ves[vv]))#, propInTrial=pits[pp], immunoDelay=immDs[dd]))
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
##        }
##     }
## }

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
