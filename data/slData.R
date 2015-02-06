if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
# slData.R
setwd('data')
require(gdata); require(data.table)

# Data from: https://data.hdx.rwlabs.org/dataset/rowca-ebola-cases

if(file.exists("allEbolaData.Rdata")){
	load("allEbolaData.Rdata")
}else{
	all <- read.xls("Data Ebola (Public).xlsx",sheet = 2) # takes a while
	save(all,file="allEbolaData.Rdata")	
}
all <- data.table(all)
all <- all[,-which(names(all) %in% c('Sources','Link')), with=F]
summary(all)

all
names(all)

sl <- all[Country=="Sierra Leone"]
sl[, Localite := gsub(" ", "", Localite)]
sl[, unique(Localite)]
sl[Localite=='Westernarearural', Localite := 'WesternAreaRural']
sl[Localite=='Westernareaurban', Localite := 'WesternAreaUrban']
sl[, unique(Localite)]
sl <- sl[Localite!='WesternArea']
sl <- sl[!grepl('National', Localite)]
sl[, table(Localite)]
sl[Localite=='Port', Localite := 'PortLoko']
sl[, table(Localite)]
sl <- sl[Localite!='Freetown']
sl[, table(Localite)]
sl[,Date:=as.Date(Date)]
sl[,Value:=as.numeric(levels(Value)[Value])]
sl[,Localite:=factor(Localite)]
   
sl[,unique(Category)]
sl[Category=='New cases']

pdf('../Figures/SL subnational data.pdf',  w = 12, h = 7)
par(lwd=1, mfrow = c(2,3))
shortRange <- range(sl[,Date])
threshDate <- as.Date('2014-12-10')
## shortRange <- c(as.Date('2014-12-01'),Sys.Date())
for(ii in 1:6) {
    plot(0,0, type = 'n', xlim = shortRange, ylim = c(0, max(sl[Category==levels(Category)[ii],Value],na.rm=T)), xaxt='n', xlab = '', bty = 'n', las = 1, ylab = '',
         main=  paste('Sierra Leone', levels(sl[,Category])[ii]))
    axis.Date(1,sl[,Date], las = 2, format = '%b-%d')
    sl[Category==levels(Category)[ii], lines(Date, Value, col = as.numeric(Localite)), Localite]
    if(ii== 5) legend('topleft', leg = sl[,unique(Localite)], col = as.numeric(sl[,unique(Localite)]), lty = 1, bty = 'n')
}
graphics.off()

sl <- sl[Category == 'New cases', list(Country, Localite, Value, Date)]
setnames(sl, 'Localite', 'reg')
setnames(sl, 'Value', 'cases')
df <- sl[, list(dd = diff(Date)), reg] ## mostly daily incidence
sum(df[,dd]!=1) ## still some problems, convert to daily incidence
sl[,int := c(NA, diff(Date)), reg]
sl[,inc := cases/int]

## 2014-09-16 repeatd twice, pick the version with more cases
sl[Date==sl[int==0,unique(Date)],]
sl[,dupl:=F]
sl[Date==sl[int==0,unique(Date)], dupl:=duplicated(Date), reg]
sl[,sum(dupl)]
sl <- sl[dupl==F]


pdf('../Figures/SL cleaned subnational data.pdf',  w = 6, h = 5)
xRange <- c(as.Date('2014-07-01'), as.Date('2015-02-01'))
par(lwd=1, 'ps' = 12)
ymax <- sl[,max(inc,na.rm=T)]*1.4
plot(0,0, type = 'n', xlim = xRange, ylim = c(0, ymax), xaxt='n', xlab = '', bty = 'n', las = 1, ylab = '', 
     main= 'Sierra Leone Incident Cases')
axis.Date(1, at=seq.Date(xRange[1], xRange[2], by = 'month'), las = 2, format = '%b-%d')
sl[, lines(Date, inc, col = as.numeric(reg)), reg]
legend('top', leg = sl[,unique(reg)], col = as.numeric(sl[,unique(reg)]), lty = 1, bty = 'n', ncol=3, cex = .75)
graphics.off()

save(sl, file='cleanSLData.Rdata')

