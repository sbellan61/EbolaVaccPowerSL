if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)
args=(commandArgs(TRUE)) ## load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    print(paste0('loading in ', args, ' from R CMD BATCH'))
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }}else{ ## otherwise set to some values for testing
        simNum <- 99
        nsims <- 3
        trial <- 'SWCT'
        batchdirnm <- 'Results'
    }

system.time(
    stps <- list(res = simNtrials(parms=makeParms(trial), N = nsims, verbose=0),
                 simNum=simNum, nsims=nsims, trial=trial)
    )

save(stps, file = file.path(batchdirnm, paste0('stps',simNum,'.Rdata')))
