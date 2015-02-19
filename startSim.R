if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R','ExpFit.R'), source)

args <- (commandArgs(TRUE)) ## load arguments from R CMD BATCH
print(args)
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    print(paste0('loading in ', args, ' from R CMD BATCH'))
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}else{
    parmArgs <- list(trial='RCT', vaccEff = 0, weeklyDecay = 1, weeklyDecayVar = 0, varClus = 0, sdLogIndiv=0, small=F)
    nsims <- 1
    seed <- 1
    simNum <- 99999
    saveNm <- 'simFPX-'
}

## seed=25;trial="SWCT";ord="none";propInTrial=0.03;sdLogIndiv=1;delayUnit=7;vaccEff=0;simNum=25;batchdirnm="BigResults/SLSimsSW";saveNm="simSL-3-";nsims=30;reordLag=14;nboot=20;
parmArgs <- subsArgs(as.list(environment()), makeParms)
print(parmArgs)
parms <- do.call(makeParms, parmArgs)
saveFl <- file.path(batchdirnm, paste0(saveNm,formatC(simNum, width=6, flag="0"),'.Rdata'))

system.time(sim <- simNtrials(seed=seed, parms=parms, N=nsims, flnm=saveFl, verbose = 1, verbFreq=10))
sim <- list(sim=sim, parms=parms, seed=seed, simNum=simNum)
save(sim, file = saveFl)

rm(list=ls(all=T))
gc()
