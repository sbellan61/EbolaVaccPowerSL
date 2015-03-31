####################################################################################################
## Initialize simulation runs with parameters fed in from R CMD BATCH.
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
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','ExpFit.R'), source)

args <- (commandArgs(TRUE)) ## load arguments from R CMD BATCH
print(args)
if(length(args)>0)  { ## Then cycle through each element of the list and evaluate the expressions.
    print(paste0('loading in ', args, ' from R CMD BATCH'))
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}else{ ## if not fed in command line parameters
seed=120;trial="SWCT";ord="none";propInTrial=0.3;sdLogIndiv=1;delayUnit=7;immunoDelay=21;vaccEff=0.9;remStartFin="TRUE";remProtDel="TRUE";simNum=2880;batchdirnm="BigResults/SLSimsFinalPTCorr";saveNm="simSL-bigPit-";nsims=1;reordLag=14;nboot=20;trialStartDate="2015-02-18"
}

verbose <- 1 ## increasing #'s allow browsing through all simulatin code (>50 is line by line)
parmArgs <- subsArgs(as.list(environment()), makeParms)
print(parmArgs)
parms <- do.call(makeParms, parmArgs)
saveFl <- file.path(batchdirnm, paste0(saveNm, sprintf("%06d", simNum),'.Rdata'))
modsToDo <- list('CoxME','GLMclus','GLMFclus') ##,'GLMMclusBy') #'GLMMclusFr')

system.time(sim <- simNtrials(seed=seed, parms=parms, N=nsims, flnm=saveFl, verbFreq=10))
sim <- list(sim=sim, parms=parms, seed=seed, simNum=simNum)
save(sim, file = saveFl)

rm(list=ls(all=T))
gc()

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
