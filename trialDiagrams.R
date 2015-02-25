require(ggplot2); require(data.table);

cats <- 20
eclipseT <- 3
addT <- 3

weeks <- cats + 2 + eclipseT + addT

spacing <- 0.02

initHaz <- (runif(cats)+0.5)*3
declRate <- 0.03
steadyDecline <- c(1,cumprod(rep(1-declRate, weeks-1)))
haz <- outer(steadyDecline, initHaz)
dim(haz) <- NULL

dat <- data.table(
  week=rep(seq(1,weeks), times = cats),
  cluster_id = rep(seq(1,cats), each = weeks),
  hazHomogeneous = haz
)


dat[ week > (cluster_id+1), status := "eclipse"]
dat[ week > (cluster_id+1+eclipseT), status := "vaccinated"]
dat[is.na(status), status := "unvaccinated"]

dat$status <- factor(dat$status, levels=c("unvaccinated","eclipse","vaccinated"))

p <- ggplot(dat) + theme_bw() +
  aes(xmin = week-1+spacing, xmax = week-spacing, 
      ymin = cluster_id-1+spacing, ymax = cluster_id - spacing, fill = hazHomogeneous, alpha=status) +
  scale_alpha_discrete(range=c(1,.5,.3), labels=c("unvaccinated","eclipse","vaccinated")) + scale_fill_continuous(low="blue", high="red")
p + geom_rect() ## comparison plot

dat.poly <- data.table(
  week = rep(dat$week, 3),
  cluster_id = rep(dat$cluster_id, 3),
  x = c(dat$week-1, dat$week, dat$week),
  y=c(dat$cluster_id, dat$cluster_id, dat$cluster_id+1),
  hazHomogeneous = rep(dat$hazHomogeneous,3),
  status = rep(dat$status, 3)
)

p + geom_polygon(data=dat.poly, mapping = aes(group=interaction(cluster_id, week), x=x, y=y))
## now need to get random sample of initHaz and decline rates

