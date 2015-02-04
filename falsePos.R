if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

nullParms <- makeParms(vaccEff = 0, weeklyDecay = 1, weeklyDecayVar = 0)
trialTypes <- c('RCT','FRCT','SWCT','CRCT')
nsims <- 2
seed <- 1

for(tri in trialTypes)
    assign(paste0('sim',tri), simNtrials(seed = seed, parms = within(nullParms, {trial=tri}) , N=nsims, verbose = 0))

sims <- NULL
for(tri in trialTypes) sims <- rbind(sims, c(trial = tri, get(paste0('sim',tri))))
as.matrix(sims)

## system.time(print(sim <- simNtrials(parms=makeParms('RCT', clusSize=300, vaccEff = .8), N=5, verbose = 0)))

