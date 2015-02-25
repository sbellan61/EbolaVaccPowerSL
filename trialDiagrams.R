require(ggplot2); require(data.table);

cats <- 20
eclipseT <- 3
addT <- 3

weeks <- cats + 2 + eclipseT + addT

weekFact <- factor(seq(1,weeks))
catFact <- factor(seq(1,cats))

spacing <- 0.02

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


dat[ as.numeric(week) > (as.numeric(cluster_id)+1), status := "eclipse"]
dat[ as.numeric(week) > (as.numeric(cluster_id)+1+eclipseT), status := "vaccinated"]
dat[is.na(status), status := "unvaccinated"]

dat$status <- factor(dat$status, levels=c("eclipse","vaccinated", "unvaccinated"))

p <- ggplot(dat) +
  aes(x=week, y=cluster_id, xmin = as.numeric(week)-0.5+spacing, xmax = as.numeric(week)+0.5-spacing, 
      ymin = as.numeric(cluster_id)-0.5+spacing, ymax = as.numeric(cluster_id)+0.5 - spacing, fill = hazHeterogeneous) +
  scale_fill_continuous(low="blue", high="red") +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + labs(y="cluster id", fill="hazard")
p + scale_alpha_manual(values=c(0.7, 0.4, 1)) + geom_rect(aes(alpha=status, fill=hazHomogeneous)) +
  labs(title="SWT, homogeneous hazard shapes")  ## comparison plot

p + scale_alpha_manual(values=c(0.7, 0.4)) +
  geom_rect(data=dat[status == "unvaccinated"]) +
  geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymax=as.numeric(cluster_id), alpha=status)) +
  geom_rect(data=dat[status != "unvaccinated"], mapping = aes(ymin=as.numeric(cluster_id))) +
  labs(title="RCT, inhomogeneous hazard shapes")
## now need to get random sample of initHaz and decline rates

