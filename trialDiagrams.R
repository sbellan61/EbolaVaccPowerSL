require(ggplot2); require(data.table); require(gridExtra)
source('multiplot.R')

set.seed(6)

cats <- 20#length(unique(ht$cluster))
eclipseT <- 3
addT <- 3

## weeks <- cats + 2 + eclipseT + addT
weeks <- 24

weekFact <- seq(-1, weeks)
catFact <- factor(LETTERS[seq(1,cats)])

yspacing <- 0.02
xspacing <- 0

initHaz <- (runif(cats,0,4)+0.5)*4
declRate <- 0.06
steadyDecline <- c(1,cumprod(rep(1-declRate, weeks+1)))
varyingDecline <- function() c(1, cumprod(rep(1-runif(1,max=5)*declRate, weeks + 1)))
hazHo <- outer(steadyDecline, initHaz)
dim(hazHo) <- NULL
hazHe <- sapply(initHaz, function(ih) ih*varyingDecline())
dim(hazHe) <- NULL

states <- c("unvaccinated", "protective delay","vaccinated")

dat <- data.table(
  week=rep(weekFact, times = cats),
  cluster_id = rep(catFact, each = weeks+2),
  hazHomogeneous = hazHo,
  hazHeterogeneous = hazHe,
  status = factor("unvaccinated", levels=states),
  order_status = factor("unvaccinated", levels=states),
  frct_order_status = factor("unvaccinated", levels=states),
  frct_status = factor("unvaccinated", levels=states)
)
setkey(dat, cluster_id, week)

stateFact <- factor(states, levels=states)

dat[ week >= as.numeric(cluster_id), status := "protective delay"]
dat[ week >= as.numeric(cluster_id)/2, frct_status := "protective delay"]

dat[ week >= as.numeric(cluster_id)+eclipseT, status := "vaccinated"]
dat[ week >= as.numeric(cluster_id)/2+eclipseT, frct_status := "vaccinated"]

for (w in 1:cats) {
  excludeclusters <- dat[(week == w) & (order_status != "unvaccinated"), ]$cluster_id
  markcluster <- dat[(week == (w-2)) & !(cluster_id %in% excludeclusters) & (hazHeterogeneous == dat[(week == (w-2)) & !(cluster_id %in% excludeclusters), max(hazHeterogeneous)]),]$cluster_id
  dat[(week >= w) & cluster_id == markcluster, order_status := "protective delay"]
  dat[(week >= (w+eclipseT)) & cluster_id == markcluster, order_status := "vaccinated"]
}

frct_clusters <- 2

for (w in 1:round(cats/2)) {
  for (i in 1:frct_clusters) {
    excludeclusters <- dat[(week == w) & (frct_order_status != "unvaccinated"), ]$cluster_id
    markcluster <- dat[(week == (w-2)) & !(cluster_id %in% excludeclusters) & (hazHeterogeneous == dat[(week == (w-2)) & !(cluster_id %in% excludeclusters), max(hazHeterogeneous)]),]$cluster_id
    dat[(week >= w) & cluster_id == markcluster, frct_order_status := "protective delay"]
    dat[(week >= (w+eclipseT)) & cluster_id == markcluster, frct_order_status := "vaccinated"]
  }
}

frctvaxord <- dat[frct_order_status != "unvaccinated", list(vaxorder = min(week)), by="cluster_id"]
frctvaxord[,week := vaxorder]
setkey(frctvaxord, cluster_id, week)
joined <- dat[frctvaxord, list(hazHeterogeneous, vaxorder, cluster_id)]
setkey(joined, vaxorder, hazHeterogeneous)
joined[,frct_vaxorder := .I]

#frctvaxord[dat, hazHeterogeneous, by=c("vaxorder")]

vaxord <- dat[order_status != "unvaccinated", list(vaxorder = min(week)), by="cluster_id"]
setkey(vaxord, "vaxorder")
vaxord$vaxorder <- factor(vaxord$vaxorder-min(vaxord$vaxorder)+1)

dat <- merge(dat, vaxord, by="cluster_id")
dat <- merge(dat, joined[,frct_vaxorder,by=cluster_id], by="cluster_id")
dat <- dat[week > 0,]

p <- ggplot(dat) +
  aes(x=week, y=cluster_id, xmin = week-0.5+xspacing, xmax = week+0.5-xspacing, 
      ymin = as.numeric(cluster_id)-0.5+yspacing, ymax = as.numeric(cluster_id)+0.5 - yspacing, fill = hazHeterogeneous) +
  scale_x_discrete(breaks=1:weeks, labels={ temp <- 1:weeks; temp[temp %% 4 != 0] <- ''; temp }, name="") +
  scale_alpha_manual(guide=F, values = c(1, 0.5, 0.4), drop = F, name="participant status") +
      scale_y_discrete(labels ='', name="") +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
  )  + labs(y="cluster", fill="hazard") 

pnleg <-   p + scale_fill_continuous(guide = F, low="dodger blue", high=rgb(.5,0,0), breaks=function(lims) {
    lim <- as.numeric(lims)
    del <- 0.05*(lim[2]-lim[1])
    c(lim[1]+del, lim[2]-del)
}, labels=c("low","high"), name="infection hazard") 

pleg <-   p + scale_fill_continuous(low="dodger blue", high=rgb(.5,0,0), breaks=function(lims) {
    lim <- as.numeric(lims)
    del <- 0.05*(lim[2]-lim[1])
    c(lim[1]+del, lim[2]-del)
}, labels=c("low","high"), name="infection hazard") 


OrigOrderRCTtu <- pnleg +  
  geom_rect(data=dat[order_status == "unvaccinated"]) +
  geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=order_status)) +
  geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
  labs(title="risk-prioritized RCT")
OrigOrderRCTtu

SWCT <- pnleg +  geom_rect(aes(alpha=status)) + labs(alpha="") +
    labs(title="SWCT")
SWCT

SWCTtu <- pnleg + aes(y=vaxorder, ymin = as.numeric(vaxorder)-0.5+yspacing, ymax = as.numeric(vaxorder)+0.5 - yspacing) +
  geom_rect(aes(alpha=status)) + labs(alpha="") +
  labs(title="SWCT, underlying order by risk priority")
SWCTtu

RCTtu <- pnleg +  aes(y=vaxorder, ymin = as.numeric(vaxorder)-0.5+yspacing, ymax = as.numeric(vaxorder)+0.5 - yspacing) +
    geom_rect(data=dat[order_status == "unvaccinated"]) +
    geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymax=as.numeric(vaxorder), alpha=order_status)) +
    geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymin=as.numeric(vaxorder))) +
    labs(title="risk-prioritized RCT") +  labs(alpha="")
RCTtu

SimInst <- pnleg +  aes(y=vaxorder, ymin = as.numeric(vaxorder)-0.5+yspacing, ymax = as.numeric(vaxorder)+0.5 - yspacing) +
  geom_rect(mapping = aes(ymax=as.numeric(vaxorder), alpha=factor(ifelse(week <= 3,"protective delay","vaccinated"), levels=states) )) +
  geom_rect(mapping = aes(ymin=as.numeric(vaxorder))) +
  labs(title="Simultaneous Instant RCT, by Hazard") +  labs(alpha="")
SimInst

SimInstByID <- pnleg +
  geom_rect(mapping = aes(ymax=as.numeric(cluster_id), alpha=factor(ifelse(week <= 3,"protective delay","vaccinated"), levels=states) )) +
  geom_rect(mapping = aes(ymin=as.numeric(cluster_id))) +
  labs(title="Simultaneous Instant RCT, by Cluster ID") +  labs(alpha="")
SimInstByID

SWCThom <- pnleg +  geom_rect(aes(alpha=status, fill=hazHomogeneous)) +
    labs(title="SWCT, homogeneous hazard shapes")
SWCThom

RCTnone <- pnleg +  
    geom_rect(data=dat[status == "unvaccinated"]) +
    geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
    geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=status)) +
    labs(title="random ordered RCT")
RCTnone

FRCT <- pnleg +
  geom_rect(data=dat[frct_status == "unvaccinated"]) +
  geom_rect(data=dat[frct_status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
  geom_rect(data=dat[frct_status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=frct_status)) +
  labs(title="random ordered FRCT")
FRCT

FRCTByHaz <- pnleg + aes(y=vaxorder, ymin = as.numeric(vaxorder)-0.5+yspacing, ymax = as.numeric(vaxorder)+0.5 - yspacing) +
  geom_rect(data=dat[frct_status == "unvaccinated"]) +
  geom_rect(data=dat[frct_status != "unvaccinated"], mapping = aes(ymin=as.numeric(vaxorder))) +
  geom_rect(data=dat[frct_status != "unvaccinated"], mapping = aes(ymax=as.numeric(vaxorder), alpha=frct_status)) +
  labs(title="random ordered FRCT, underlying order by hazard")
FRCTByHaz

FRCTHazOrdByHaz <- pnleg + aes(y=frct_vaxorder, ymin = as.numeric(frct_vaxorder)-0.5+yspacing, ymax = as.numeric(frct_vaxorder)+0.5 - yspacing) +
  geom_rect(data=dat[frct_order_status == "unvaccinated"]) +
  geom_rect(data=dat[frct_order_status != "unvaccinated"], mapping = aes(ymin=as.numeric(frct_vaxorder))) +
  geom_rect(data=dat[frct_order_status != "unvaccinated"], mapping = aes(ymax=as.numeric(frct_vaxorder), alpha=frct_order_status)) +
  labs(title="haz ordered FRCT, underlying order by hazard")
FRCTHazOrdByHaz

FRCTHazOrdByID <- pnleg + 
  geom_rect(data=dat[frct_order_status == "unvaccinated"]) +
  geom_rect(data=dat[frct_order_status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
  geom_rect(data=dat[frct_order_status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=frct_order_status)) +
  labs(title="haz ordered FRCT, underlying order by ID")
FRCTHazOrdByID

pdf('Figures/Fig 3 schematic.pdf', w = 6.5, h = 4)
multiplot(SWCT, OrigOrderRCTtu, cols = 2)
graphics.off()


pdf('Figures/Fig SX schematic.pdf', w = 6.5, h = 7)
multiplot(SWCT, RCTnone, OrigOrderRCTtu, FRCTHazOrdByID, cols = 2)
graphics.off()
