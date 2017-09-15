#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation- 8 Day Chlorophyll")
df<-read.csv("subsetTracks30MinInterpolationLandMod.csv",header=T,sep=",")

# rename columns
names(df)[names(df) == 'MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A'] <- 'chloro'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'
names(df)[names(df) == 'tag.local.identifier'] <- 'ID'

# remove row with missing chlorophyll value 
df<-df[complete.cases(df[ , 11]),]

# convert times from factors to POSIX
df$date<-as.POSIXct(df$date, format= "%d/%m/%Y %H:%M", tz = "GMT")
head(df)

# split up data by location - Ireland or Scotland 
df <- df[df$location=="Scotland" , ] 
df<-droplevels(df)

# prepare data with moveHMM
stormData <- df[,c(4,5,8,11,13)]
head(stormData)
stormData <- prepData(stormData,type="LL",coordNames=c("lon","lat"))
# plot(stormData,compact=T)

# add cosinor covariate based on hour of day
stormData$hour <- as.integer(strftime(stormData$date, format = "%H", tz="GMT"))
head(stormData)

# label 2 states
stateNames <- c("exploratory", "encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(10,5,1,2),angle=c(10,5)) # it's mean1,mean2,sd1,sd2 for step lengths

# fit model
m1 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m1

# formula for transition probabilities
formula <- ~ chloro * cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m2

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ chloro * cosinor(hour, period = 24),
                       sd = ~ chloro * cosinor(hour, period = 24)),
           angle = list(concentration = ~ chloro))
# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)

# fit model
m3 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
             formula = formula)

# label 3 states
stateNames <- c("exploratory", "encamped", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m4 <- list(step=c(10,5,1,1,2,1),angle=c(10,5,1))

# fit model
m4 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m4

# formula for transition probabilities
formula <- ~ chloro * cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m5 <- getPar0(model=m4, formula=formula)
# fit model
m5 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m5$Par,
             beta0=Par0_m5$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m5

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ chloro * cosinor(hour, period = 24),
                       sd = ~ chloro * cosinor(hour, period = 24)),
           angle = list(concentration = ~ chloro))
# initial parameters (obtained from nested model m5)
Par0_m6 <- getPar0(model=m5, formula=formula, DM=DM)

# fit model
m6 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m6$Par,
             beta0 = Par0_m6$beta, DM = DM, stateNames = stateNames,
             formula = formula)

plot(m6,covs = data.frame(hour=6))

# compute the pseudo-residuals
pr <- pseudoRes(m6)
# time series, qq-plots, and ACF of the pseudo-residuals
plotPR(m6,title("m6"))

AIC(m1,m2,m3,m4,m5,m6)

