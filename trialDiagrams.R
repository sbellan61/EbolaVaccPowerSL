require(ggplot2); require(data.table);

cats <- 20
eclipseT <- 3
addT <- 3

weeks <- cats + 2 + eclipseT + addT

weekFact <- factor(seq(1,weeks), levels=seq(1,weeks))
catFact <- factor(seq(1,cats), levels=seq(1,cats))

yspacing <- 0.02
xspacing <- 0

initHaz <- (runif(cats)+0.5)*4
declRate <- 0.03
steadyDecline <- c(1,cumprod(rep(1-declRate, weeks-1)))
varyingDecline <- function() c(1, cumprod(rep(1-runif(1,max=5)*declRate, weeks - 1)))
hazHo <- outer(steadyDecline, initHaz)
dim(hazHo) <- NULL
hazHe <- sapply(initHaz, function(ih) ih*varyingDecline())
dim(hazHe) <- NULL

dat <- data.table(
  week=rep(weekFact, times = cats),
  cluster_id = rep(catFact, each = weeks),
  hazHomogeneous = hazHo,
  hazHeterogeneous = hazHe
)

states <- c("unvaccinated", "seroconversion lag","vaccinated")

stateFact <- factor(states, levels=states)

dat[ as.numeric(week) > (as.numeric(cluster_id)+1), status := factor("seroconversion lag", levels=states)]
dat[ as.numeric(week) > (as.numeric(cluster_id)+1+eclipseT), status := factor("vaccinated", levels=states)]
dat[is.na(status), status := factor("unvaccinated", levels=states)]

dat$status <- factor(dat$status, levels=states)

dat$order_status <- factor("unvaccinated", levels=states)

for (w in (1:cats)+2) {
  excludeclusters <- dat[(week == w) & (order_status != "unvaccinated"), ]$cluster_id
  markcluster <- dat[(week == (w-2)) & !(cluster_id %in% excludeclusters) & (hazHeterogeneous == dat[(week == (w-2)) & !(cluster_id %in% excludeclusters), max(hazHeterogeneous)]),]$cluster_id
  dat[(as.numeric(week) >= w) & cluster_id == markcluster, order_status := factor("seroconversion lag", levels=states)]
  dat[(as.numeric(week) >= (w+eclipseT)) & cluster_id == markcluster, order_status := factor("vaccinated", levels=states)]
}

vaxord <- dat[order_status != "unvaccinated", list(vaxorder = min(as.numeric(week))), by="cluster_id"]
setkey(vaxord, "vaxorder")
vaxord$vaxorder <- factor(vaxord$vaxorder-min(vaxord$vaxorder)+1)
dat <- merge(dat, vaxord, by="cluster_id")

p <- ggplot(dat) +
  aes(x=week, y=cluster_id, xmin = as.numeric(week)-0.5+xspacing, xmax = as.numeric(week)+0.5-xspacing, 
      ymin = as.numeric(cluster_id)-0.5+yspacing, ymax = as.numeric(cluster_id)+0.5 - yspacing, fill = hazHeterogeneous) +
  scale_fill_continuous(low="blue", high="red", breaks=function(lims) {
    lim <- as.numeric(lims)
    del <- 0.05*(lim[2]-lim[1])
    c(lim[1]+del, lim[2]-del)
  }, labels=c("low","high")) +
  scale_x_discrete(limits=weekFact) +
  scale_alpha_manual(values = c(1, 0.7, 0.4), drop = F) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + labs(y="cluster id", fill="hazard")

p + geom_rect(aes(alpha=status, fill=hazHomogeneous)) +
  labs(title="SWT, homogeneous hazard shapes")

p + geom_rect(aes(alpha=status)) +
  labs(title="SWT, heterogeneous hazard shapes")

p + 
  geom_rect(data=dat[status == "unvaccinated"]) +
  geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
  geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=status)) +
  labs(title="RCT, heterogeneous hazard shapes")

p + aes(y=vaxorder, ymin = as.numeric(vaxorder)-0.5+yspacing, ymax = as.numeric(vaxorder)+0.5 - yspacing) +
  scale_y_discrete(labels = vaxord$cluster_id, name="cluster id") +
  geom_rect(data=dat[order_status == "unvaccinated"]) +
  geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymax=as.numeric(vaxorder), alpha=order_status)) +
  geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymin=as.numeric(vaxorder))) +
  labs(title="RCT, heterogeneous hazard shapes, clusters ordered by hazard", alpha="status")

## use hazard-based approach to vax clusters, re-order according to week of vaccination

## now need to get random sample of initHaz and decline rates

