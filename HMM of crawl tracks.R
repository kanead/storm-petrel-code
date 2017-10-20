#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
library(momentuHMM)
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\crawl SP-Movebank")
stormData <- read.csv("crawl SP-3493787305078051578.csv",header = T , sep = ",")
head(stormData)

stormData<-stormData[,c("location.long","location.lat","timestamp", "individual.local.identifier","ETOPO1.Elevation",
              "NASA.Distance.to.Coast", "MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A")]
names(stormData)[names(stormData) == 'location.long'] <- 'lon'
names(stormData)[names(stormData) == 'location.lat'] <- 'lat'
names(stormData)[names(stormData) == 'individual.local.identifier'] <- 'ID'
names(stormData)[names(stormData) == 'ETOPO1.Elevation'] <- 'bath'
names(stormData)[names(stormData) == 'NASA.Distance.to.Coast'] <- 'coast'
names(stormData)[names(stormData) == 'MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A'] <- 'chloro'
head(stormData)

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(stormData[,1:2],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm + ellps=WGS84")) # 29 = IRE or 30 = UK
# add UTM locations to data frame
stormData$x <- attr(utmcoord,"coords")[,1]
stormData$y <- attr(utmcoord,"coords")[,2]

plot(stormData$x , stormData$y)

stormData <- momentuHMM::prepData(stormData)
###################################################################################
# label states
stateNames <- c("exploratory","encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))

# fit model
m1 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 5)

m1
###################################################################################
# formula for 2 state model with transition probabilities 
formula <- ~ bath 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m2
###################################################################################
# label 3 states
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))

# fit model
m3 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m3
###################################################################################
# formula for 3 state model with transition probabilities
formula <- ~ bath 
# initial parameters (obtained from nested model m1)
Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
m4 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
             beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1)

###################################################################################
AIC(m1,m2,m3,m4)
###################################################################################
# Chlorophyll covariate 
stormData<-head(stormData,-5)
###################################################################################
# formula for 2 state model with transition probabilities 
stateNames <- c("exploratory","encamped")
formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m5 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
                         beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m5
###################################################################################
# formula for 3 state model with transition probabilities
stateNames <- c("transiting", "foraging", "resting")
formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
m6 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
                         beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m6
###################################################################################
AIC(m1,m2,m3,m4,m5,m6)
###################################################################################
# Chlorophyll covariate 
stormData<-head(stormData,-5)
###################################################################################
# formula for 2 state model with transition probabilities 
stateNames <- c("exploratory","encamped")
formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m5 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
                         beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m5
###################################################################################
# formula for 3 state model with transition probabilities
stateNames <- c("transiting", "foraging", "resting")
formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
m6 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
                         beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m6
###################################################################################
AIC(m1,m2,m3,m4,m5,m6,m7,m8)
###################################################################################