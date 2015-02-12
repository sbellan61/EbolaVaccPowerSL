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
    parmArgs <- subsArgs(as.list(environment()), makeParms)
}else{
    parmArgs <- list(trial='RCT', vaccEff = 0, weeklyDecay = 1, weeklyDecayVar = 0, varClus = 0, sdLogIndiv=0, small=F)
    nsims <- 1
    seed <- 1
    simNum <- 99999
    saveNm <- 'simFPX-'
}

print(parmArgs)
parms <- do.call(makeParms, parmArgs)

system.time(sim <- simNtrials(seed=seed, parms=parms, N=nsims, verbose = 1, verbFreq=10))
sim <- list(sim=sim, parms=parms, seed=seed)
save(sim, file = file.path(batchdirnm, paste0(saveNm,formatC(simNum, width=6, flag="0"),'.Rdata')))

rm(list=ls(all=T))
gc()
## for(tri in trialTypes)
##     assign(paste0('sim',tri), simNtrials(seed = seed, parms = within(nullParms, {trial=tri}),
##                                          N=nsims, verbose = 0, returnAll=T, showSeqStops=T, flnm = paste0('sStop',tri)))

## resfull <- seqStop(res, fullSeq = T)
## showSeqStop(resfull, paste('test',ii), width = 8, height = 6)
## system.time(print(sim <- simNtrials(parms=makeParms('RCT', clusSize=300, vaccEff = .8), N=5, verbose = 0)))

