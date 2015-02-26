require(ggplot2); require(data.table); require(gridExtra)
source('multiplot.R')

cats <- 20
eclipseT <- 3
addT <- 3

## weeks <- cats + 2 + eclipseT + addT
weeks <- 24

weekFact <- factor(seq(1,weeks))
catFact <- factor(LETTERS[seq(1,cats)])

yspacing <- 0.02
xspacing <- 0

initHaz <- (runif(cats,0,4)+0.5)*4
declRate <- 0.05
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
  hazHeterogeneous = hazHe,
  status = factor("unvaccinated", levels=states),
  order_status = factor("unvaccinated", levels=states)
)

states <- c("unvaccinated", "protective delay","vaccinated")

stateFact <- factor(states, levels=states)

dat[ as.numeric(week) > (as.numeric(cluster_id)+1), status := "protective delay"]
dat[ as.numeric(week) > (as.numeric(cluster_id)+1+eclipseT), status := "vaccinated"]

for (w in (1:cats)+2) {
  excludeclusters <- dat[(week == w) & (order_status != "unvaccinated"), ]$cluster_id
  markcluster <- dat[(week == (w-2)) & !(cluster_id %in% excludeclusters) & (hazHeterogeneous == dat[(week == (w-2)) & !(cluster_id %in% excludeclusters), max(hazHeterogeneous)]),]$cluster_id
  dat[(as.numeric(week) >= w) & cluster_id == markcluster, order_status := "protective delay"]
  dat[(as.numeric(week) >= (w+eclipseT)) & cluster_id == markcluster, order_status := "vaccinated"]
}

vaxord <- dat[order_status != "unvaccinated", list(vaxorder = min(as.numeric(week))), by="cluster_id"]
setkey(vaxord, "vaxorder")
vaxord$vaxorder <- factor(vaxord$vaxorder-min(vaxord$vaxorder)+1)
dat <- merge(dat, vaxord, by="cluster_id")

labs <- 0:24
labs[labs %% 4 !=0] <- ''

p <- ggplot(dat) +
  aes(x=week, y=cluster_id, xmin = as.numeric(week)-0.5+xspacing, xmax = as.numeric(week)+0.5-xspacing, 
      ymin = as.numeric(cluster_id)-0.5+yspacing, ymax = as.numeric(cluster_id)+0.5 - yspacing, fill = hazHeterogeneous) +
  scale_fill_continuous(low="dodger blue", high="red", breaks=function(lims) {
    lim <- as.numeric(lims)
    del <- 0.05*(lim[2]-lim[1])
    c(lim[1]+del, lim[2]-del)
  }, labels=c("low","high")) +
  scale_x_discrete(limits=weekFact, breaks = 0:24, labels = labs) +
  scale_alpha_manual(values = c(1, 0.7, 0.4), drop = F) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + labs(y="cluster", fill="hazard") 



SWT <- p + geom_rect(aes(alpha=status)) + labs(alpha="") +
    labs(title="SWT")
SWT

RCTtu <- p + aes(y=vaxorder, ymin = as.numeric(vaxorder)-0.5+yspacing, ymax = as.numeric(vaxorder)+0.5 - yspacing) +
    scale_y_discrete(labels = vaxord$cluster_id, name="cluster") +
    geom_rect(data=dat[order_status == "unvaccinated"]) +
    geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymax=as.numeric(vaxorder), alpha=order_status)) +
    geom_rect(data=dat[order_status != "unvaccinated"], mapping = aes(ymin=as.numeric(vaxorder))) +
    labs(title="RCT") +  labs(alpha="")


SWThom <- p + geom_rect(aes(alpha=status, fill=hazHomogeneous)) +
    labs(title="SWT, homogeneous hazard shapes")

RCTnone <- p + 
    geom_rect(data=dat[status == "unvaccinated"]) +
    geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
    geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=status)) +
    labs(title="RCT, heterogeneous hazard shapes")

##pdf('Figures/Fig 1 schematic.pdf', w = 3.5, h = 9)
png('Figures/Fig 1 schematic.png', w = 3.5, h = 9, units='in', res = 200)
grid_arrange_shared_legend(SWT, RCTtu, ncol = 1)
graphics.off()

