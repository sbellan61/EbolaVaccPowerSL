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
  cluster_id = as.factor(rep(seq(1,cats), each = weeks)),
  hazHomogeneous = haz
)
p <- ggplot(dat) + theme_bw() +
  aes(xmin = week-1+spacing, xmax = week-spacing, 
      ymin = as.numeric(cluster_id)-1+spacing, ymax = as.numeric(cluster_id) - spacing, fill=hazHomogeneous) +
  geom_rect()
