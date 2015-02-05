if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

batchdirnm <- file.path('Results','FalsePosSims')
tnms <- c('SWCT','RCT','FRCT','CRCT')
numEach <- 12
parmsMat <- data.frame(
    nsims = 200
  , seed =  rep(1:numEach, length(tnms))
  , trial = rep(tnms, each = numEach)
  , batchdirnm = batchdirnm
  , varClus = 0
  , weeklyDecay = 1
  , weeklyDecayVar = 0
  , vaccEff = 0
)
parmsMat$simNum <- 1:nrow(parmsMat)


fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

paste0("Results/FalsePosSims/simFP", 1:48,".Rdata") %in% fls

smp <- NULL
for(tri in tnms) {
    simNms <- paste0("Results/FalsePosSims/simFP", which(parmsMat$trial==tri),".Rdata")
    simNms <- fls[fls %in% simNms]
    tmp <- NULL
    for(ss in simNms) {
        load(ss)
        if(nrow(sim)!=200) smp <- c(smp,ss)
        tmp <- rbind(tmp, sim)
    }
    assign(paste0('wa',tri), as.data.table(tmp))
}
smp

tnms
for(tri in tnms) print(get(paste0('wa',tri))[, mean(stopped & vaccGood)])
for(tri in tnms) print(get(paste0('wa',tri))[, mean(stopped & !vaccGood)])
for(tri in tnms) print(get(paste0('wa',tri))[stopped & vaccGood, mean(stopDay)])
for(tri in tnms) print(get(paste0('wa',tri))[stopped & vaccGood, mean(ptRatioCVXimmGrpEnd)])
for(tri in tnms) print(get(paste0('wa',tri))[stopped & vaccGood, mean(caseCXimmGrpEnd + caseVXimmGrpEnd)])
names(waSWCT)
for(tri in tnms) print(dim(get(paste0('wa',tri))))



