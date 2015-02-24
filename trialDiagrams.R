require(ggplot2); require(data.table);

weeks <- 10
cats <- weeks - 2
spacing <- 0.02

initHaz <- (runif(cats)+0.5)*3
declRate <- 0.01
steadyDecline <- c(1,cumprod(rep(1-declRate, weeks-1)))
haz <- outer(steadyDecline, initHaz)
dim(haz) <- NULL

dat <- data.table(
  week=rep(seq(1,weeks), times = cats),
  cluster_id = rep(seq(1,cats), each = weeks),
  hazHomogeneous = haz
)

dat[ week > (cluster_id+1), status := "vaccinated"]
dat[is.na(status), status := "unvaccinated"]

p <- ggplot(dat) + theme_bw() +
  aes(xmin = week-1+spacing, xmax = week-spacing, 
      ymin = cluster_id-1+spacing, ymax = cluster_id - spacing, fill = hazHomogeneous, alpha=status) +
  geom_rect() + scale_alpha_discrete(range=c(1,.7)) + scale_fill_continuous(low="blue", high="red")
