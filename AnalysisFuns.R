## Construct survival data from waiting times
makeSurvDat <- function(parms) within(parms, {
    pop$immuneDayThink <- pop[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    ## pre-immunity table
    stPre <- copy(pop) # st = survival table
    stPre$startDay <- 0
    stPre$endDay <- stPre[, pmin(immuneDayThink, infectDay)]
    stPre$infected <- stPre[ ,as.numeric(infectDay <= immuneDayThink)]
    stPre$immuneGrp <-  0     ## immuneGrp is variable used for analysis, not omnietient knowledge of vaccination/immune status
    stPre <- stPre[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    ## post-immunity table
    stPost <- copy(pop)[infectDay > immuneDayThink,]
    stPost$startDay <- stPost[,immuneDayThink]
    stPost$endDay   <-  stPost[,infectDay]
    stPost$infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
    stPost$immuneGrp <- 1
    stPost <- stPost[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, immuneDayThink, startDay, endDay, infected, immuneGrp)]
    st <- rbind(stPre, stPost) ## combine tables
    rm(stPre, stPost)
}) ## careful with modifying parms, st depends on analysis a bit too (immunoDelayThink), so we can have different st for same popH

## Select subset of survival table to analyze
activeFXN <- function(parms) within(parms, { 
    ## for SWCT or unmatched CRCT always include all clusters in analysis because unvaccinated
    ## clusters are still considered to provide useful information from baseline
    st$firstActive <- 0
    if(!includeAllControlPT) { ## remove person-time observed prior to post-refractory period from data
        if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
            st[, firstActive := min(immuneDayThink), by = pair]
        if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
            st[, firstActive := min(immuneDayThink), by = cluster]
    }
    st <- st[!endDay <= firstActive] ## remove inactive observation intervals
    st[startDay < firstActive, startDay := firstActive] ## set accumulation of person time as when the cluster/pair is active
})
## p1 <- simTrial(makeParms('RCT', ord='BL', small=F), br=F)
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## s1$st[idByClus%in%1:2, list(indiv, cluster, pair, idByClus,immuneDayThink, startDay,endDay)]

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = 6*30) with(parms, {
    intervalNotStarted <- st[,startDay] > censorDay
    st <- st[!intervalNotStarted,] 
    noInfectionBeforeCensor <- st[,endDay] > censorDay
    st[noInfectionBeforeCensor, infected:=0]
    st[noInfectionBeforeCensor, endDay:=censorDay]
    st[,perstime := (endDay-startDay)]
    st <- st[perstime > 0,] 
    return(st)
})

summTrial <- function(st) list(summarise(group_by(st, cluster), sum(infected))
                               , summarise(group_by(st, cluster, immuneGrp), sum(infected))
                               , summarise(group_by(st, immuneGrp), sum(infected))
                               )

compileStopInfo <- function(minDay, vaccEffEst, tmp) {
    out <- c(stopDay=minDay, vaccEffEst, caseVacc = tmp[immuneGrp==1, sum(infected)], caseCont = tmp[immuneGrp==0, sum(infected)],
             ptVacc = tmp[immuneGrp==1, sum(perstime)], ptCont = tmp[immuneGrp==0, sum(perstime)] )
    out <- c(out, hazVacc = as.numeric(out['caseVacc']/out['ptVacc']/yearToDays),
             hazCont = as.numeric(out['caseCont']/out['ptCont']/yearToDays))
    out <- c(out, ptRatio = as.numeric(out['ptCont'] / out['ptVacc']))
    return(out)
}

## ## Do a binary search for the number of infections before the stopping point is reached: this is
## ## assumed to be when 95% CI of vaccine efficacy goes above 0
## firstStop <- function(parms, minDay=min(parms$pop$immuneDayThink) + 30, maxDay=365, verbose = 0) { ## using days to facilitate easier rounding
##     if(verbose>=2) browser()
##     midDay <- floor((minDay+maxDay)/2) ## floor to days
##     tmp <- censSurvDat(parms, midDay)
##     vaccEffEst <- doCoxPH(tmp) ## converting midDay to days from months
##     out <- compileStopInfo(minDay, vaccEffEst, tmp)
##     if (minDay >= maxDay) return(out)
##     pVal <- vaccEffEst['p']
##     if(verbose>0) print(signif(vaccEffEst,2))
##     goSmaller <- !is.na(pVal) & pVal<.05
##     if(goSmaller)
##         return(firstStop(parms, minDay, midDay, verbose))
##     return(firstStop(parms, midDay+1, maxDay, verbose)) ## output in months
## }

## Check whether stopping point has been reached at intervals
seqStop <- function(parms, start = parms$immunoDelayThink + 14, checkIncrement = 7, verbose = 0, maxDay = parms$maxInfectDay) {
    trialOngoing <- T
    checkDay <- start
    first <- T
    while(trialOngoing) {
        if(verbose>1) browser()
        if(verbose>0) print(checkDay)
        tmp <- censSurvDat(parms, checkDay)
        vaccEffEst <- try(doCoxPH(tmp), silent=F) ## converting midDay to days from months
        ## if cox model has enough info to converge check for stopping criteria
        if(!inherits(vaccEffEst, 'try-error') & !is.nan(vaccEffEst['p'])) { 
            newout <- compileStopInfo(checkDay, vaccEffEst, tmp) 
            if(first) out <- newout else out <- rbind(out, newout)
            first <- F
            if(!is.na(newout['p']))
                if(newout['p'] < .05)
                    trialOngoing <- F
            if(checkDay > maxDay) trialOngoing <- F
        }
        checkDay <- checkDay + 7
    }
    rownames(out) <- NULL
    out <- as.data.table(out)
    parms$weeklyAns <- out
    parms$endTrialDay <- tail(out$stopDay,1)
    parms$vaccEffEst <- vaccEffEst
    return(parms)
}

## Vaccinate control groups once vaccine efficacy identified
vaccEndTrial <- function(parms) within(parms, {
    EVpopH <- copy(popH)
    toChange <- EVpopH[, vaccDay==Inf]
    lastPreassignedVaccDay <- popH[vaccDay!=Inf, max(vaccDay)]
    startVaccContDay <- max(lastPreassignedVaccDay, endTrialDay)
    ## SWCT is fastest possible vaccination already, so no change

    if(trial %in% c('RCT','FRCT')) { ## remaining half get vaccinated in each cluster starting delayUnit day decision made to vaccinate controls
        vaccClusters <- EVpopH[vaccDay < endTrialDay + delayUnit & vaccDay!=Inf , unique(cluster)]            
        notYetVaccClusters <- EVpopH[vaccDay >= endTrialDay + delayUnit & vaccDay!=Inf , unique(cluster)]
        ## ##################################################
        ## OPTION 1: Continue with cluster sequence, but vaccinate full clusters. Then go back and vacc partially vacc clusters.
        if(F) {# RCTendOption==1) {
            ## vaccinate everyone (not just half) in those clusters that haven't yet had anyone vaccinated moving forward
            vaccDaysLeft <- daySeq[daySeq > lastPreassignedVaccDay] ## remaining vaccination days (delayUnit apart)
            EVpopH[cluster %in% notYetVaccClusters & vaccDay==Inf,
                   vaccDay := delayUnit*(cluster-1)] ## same as those pre-assigned to be vaccinated
            EVpopH[cluster %in% notYetVaccClusters, unique(vaccDay)]
            ## ##################################################
            ## Vaccinate the rest of partially vaccinated clusters afterwards. First get order in which to do it, if necessary.
            if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVclusIncRank <- EVpopH[cluster %in% vaccClusters & vaccDay == Inf & idByClus==1 & day == endTrialDay,
                                        rev(order(clusHaz))]
            if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                EVclusIncRank <- NULL
                for(ii in 1:length(vaccClusters)) {
                    dd <- vaccDaysLeft[ii]
                    tmpRank <- EVpopH[cluster %in% vaccClusters & idByClus==1 & day == dd,
                                      rev(order(clusHaz))]
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }
            ## Then update the vaccination days based on the order
            if(ord == 'none')
                EVpopH[cluster %in% vaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
            if(ord == 'BL') { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }
            if(ord == 'TU') {
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
            }
            rm(ii) }
        ## EVpopH[idByClus %in% c(1,1+clusSize/2) & day == daySeq[which.min(abs(daySeq - startVaccContDay))],
        ##        list(cluster, clusHaz, day, vaccDay)] ## check that it works
    }
    ## ##################################################
    ## Option 2. Vacc partially vacc clusters first, then vaccinate entire unvacc clusters.
    ## vaccinate the rest of partially vaccinated clusters afterwards
    if(F) {# RCTendOption==2) {
        ## Vaccinate the rest of partially vaccinated clusters first. First get order in which to do it, if necessary.
        vaccDaysLeft <- daySeq[daySeq > endTrialDay] ## remaining vaccination days (delayUnit apart)
            if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVclusIncRank <- EVpopH[cluster %in% vaccClusters & vaccDay == Inf & idByClus==1 & day == endTrialDay,
                                        rev(order(clusHaz))]
            if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                EVclusIncRank <- NULL
                for(ii in 1:length(vaccClusters)) {
                    dd <- vaccDaysLeft[ii]
                    tmpRank <- EVpopH[cluster %in% vaccClusters & idByClus==1 & day == dd,
                                      rev(order(clusHaz))]
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }
            ## Then update the vaccination days for unvaccinated individuals in partially vaccinated clusters based on these orders
            if(ord == 'none')
                EVpopH[cluster %in% vaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
            if(ord == 'BL') { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }
            if(ord == 'TU') {
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
            }
        ## EVpopH[cluster %in% vaccClusters & idByClus %in% c(1,clusSize/2+1) & day == 0, ] ## check ordering works
        ## Then go back and vaccinate unvaccinated clusters
        vaccDaysLeft <- daySeq[daySeq > EVpopH[cluster %in% vaccClusters, max(vaccDay)]] ## update this vector
        ## Again get order first, for totally unvaccinated clusters
        if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when *trial ended*
            EVclusIncRank <- EVpopH[cluster %in% notYetVaccClusters & vaccDay == Inf & idByClus==1 & day == endTrialDay,
                                    rev(order(clusHaz))]
        if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
            EVclusIncRank <- NULL
            for(ii in 1:length(notYetVaccClusters)) { ## HERE********* ## 
                dd <- vaccDaysLeft[ii]
                tmpRank <- EVpopH[cluster %in% notYetVaccClusters & idByClus==1 & day == dd,
                                  rev(order(clusHaz))]
                tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
            }
        }
        ## Then vaccinate
       if(ord == 'none')
                EVpopH[cluster %in% notYetVaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
            if(ord == 'BL') { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==notYetVaccClusters[EVclusIncRank[ii]] & vaccDay==Inf,
                           vaccDay := vaccDaysLeft[ii]]
                rm(ii) }
            if(ord == 'TU') {
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==notYetVaccClusters[EVclusIncRank[ii]] & vaccDay==Inf,
                           vaccDay := vaccDaysLeft[ii]]
            }
        }
        ## ##################################################
        ## Option 3. Continue with cluster sequence, but vaccinate full clusters. *Simultaneously* begin going
        ## back and vacc partially vacc clusters.
        if(F) {# RCTendOption==3) {

        }
    if(trial=='CRCT') { ## remaining clusters get vaccinated starting delayUnit after endTrialDay or when the last cluster 
            if(ord != 'TU') ## same for no ordring & baseline ordering, remaining clusters will be ordered aleady
                EVpopH[vaccDay==Inf, vaccDay := startVaccContDay + delayUnit*(cluster - numClus/2)] ## remaining clusters go numClus/2+1 : numClus in CRCT
            if(ord == 'TU') {


            }
    ## Reset immune indices for these individuals
    EVpopH[toChange, immuneDay := vaccDay + immunoDelay] 
    EVpopH[toChange, vacc := day >= vaccDay]
    EVpopH[toChange, immune := day >= immuneDay]

    ## ## Resample infection times for vacc
    ## EVpopH[toChange & vacc==TRUE,]
    
})

t1$endTrialDay <- 70
vaccEndTrial(t1)

## * Finish the trial (but vaccinating whole clusters), then go back and
## clean up the half-vaccinated clusters

## * Clean up the half-vaccinated clusters, and then finish with the
## unvaccinated clusters

## * Assume some sort of logistical miracle that allows you to do both.


casesInTrial <- function(parms, maxDayCaseDay = 6*30) sum(with(parms$pop, infectDay < maxDayCaseDay))

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, check=F, verbose=0) {
    set.seed(seed)
    for(ii in 1:N) {
        if(verbose>1) browser()
        res <- simTrial(parms)
        stopPoint <- tail(seqStop(res),1)
        if(ii==1) out <- stopPoint else out <- rbind(out, stopPoint)
        if(check) {
            doCoxPH(censSurvDat(res$st, stopPoint$stopDay))
            doCoxPH(censSurvDat(res$st, stopPoint$stopDay+1))
        }
    }
    rownames(out) <- NULL
    return(out)#as.data.table(out))
}

simNwrp <- function(parms=makeParms(), NperCore = 10, check=F, ncores=12) {
    out <- mclapply(1:ncores, simNtrials, N = NperCore, mc.cores = ncores)
    out <- do.call(rbind.data.frame, out)
}
