# slData.R

require(gdata)

# Data from: https://data.hdx.rwlabs.org/dataset/rowca-ebola-cases

if(!file.exists("allEbolaData.Rdata")){
	load("allEbolaData.Rdata")
}else{
	all <- read.xls("Data Ebola (Public).xlsx",sheet = 2) # takes a while
	save(all,file="allEbolaData.Rdata")	
}

head(all)
dat <- subset(all,)