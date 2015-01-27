# Make a trial population with a given number of clusters of a given size. Put the people in clusters; each cluster has a beta-distributed underlying antibody proportion, which is applied to the people.
makePop <- function(numClus, clusSize){
	pop <- data.frame(indiv=as.factor(1:(numClus*clusSize))
		, cluster=as.numeric(gl(n=numClus, k=clusSize))
	)
	return(pop)
}

# Simulate the stepped wedge trial for one person
# Currently we have no individual, cluster or time effects (other than proportion immune)
# We cut off each person's time series once they get sick (using cumsum), if they do

# drop (as an argument) is purely for debugging. In real life, we don't observe people once they've had illness, we drop them.
simulateIndiv <- function(p, numClus, drop=TRUE){

	pf <- data.frame(period=1:(numClus-1))
	pf <- with(as.list(p), within(pf, {
		vacc <- ifelse(period>=p$cluster, 1, 0)
		prob <- ifelse(vacc==1, p$postProb, p$preProb)
		dis <- rbinom(length(period), size=1, p=prob)
		predis <- c(0, cumsum(dis[-length(dis)]))
	}))

	if (drop){
		pf <- subset(pf, predis==0)
	}
	pf <- subset(pf, select=c(period, vacc, dis))

	pf <- with(as.list(p), within(pf, {
		abx=abx
		indiv=indiv
		cluster=as.factor(cluster)
	}))

	return(pf)
}

# Make a population and simulate all of the stepped-wedge intervals
simulatePop <- function(numClus, clusSize, pAbx, sAbx
	, vaccProt, abxProt, suscFun, nProb, drop=TRUE
){
	people <- makePop(numClus, clusSize, pAbx, sAbx)

	# Calculate susceptibility of each person before and after
	# vaccination
	people <- within(people, {
		preProb <- nProb*suscWrapper(
			suscFun, abx, vacc=0*abx, abxProt, vaccProt)
		postProb <- nProb*suscWrapper(
			suscFun, abx, vacc=1+0*abx, abxProt, vaccProt)
	})

	outcomes <- simulateIndiv(people[1, ], numClus, drop)
	for(p in 2:nrow(people)){
		outcomes <- rbind(outcomes, 
			simulateIndiv(people[p, ], numClus))
	}
	return(outcomes)
}

# Given a model fit, return a one-tailed P value (goes from 0 for vaccine highly protective to 1 for vaccine highly risky). This is a copy from the other project, which is bad.
oneTailP <- function(m){
	s <- summary(m)
	v <- as.list(s$coef["vacc",])
	if (v$z>0) return(1-v$P/2)
	return(v$P/2)
}

rmodP <- function(people, mform){
	mod <- glmer(mform, data=people,
		family=binomial(link="cloglog")
	)

	if(numTrials<20) print(summary(mod))
 
	return(oneTailP(mod))
}

modP <- function(people, mform){
	mod <- glm(mform, data=people,
		family=binomial(link="cloglog")
	)

	if(numTrials<20) print(summary(mod))
 
	return(oneTailP(mod))
}
