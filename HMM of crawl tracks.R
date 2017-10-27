#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
library(momentuHMM)
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\crawl SP-Movebank")
stormData <- read.csv("crawl SP 5 day chloro.csv",header = T , sep = ",")
head(stormData)

stormData<-stormData[,c("location.long","location.lat","timestamp.1", "individual.local.identifier","X5dayAverage", "ETOPO1.Elevation")]
names(stormData)[names(stormData) == 'location.long'] <- 'lon'
names(stormData)[names(stormData) == 'location.lat'] <- 'lat'
names(stormData)[names(stormData) == 'individual.local.identifier'] <- 'ID'
names(stormData)[names(stormData) == 'ETOPO1.Elevation'] <- 'bath'
#names(stormData)[names(stormData) == 'NASA.Distance.to.Coast'] <- 'coast'
names(stormData)[names(stormData) == 'X5dayAverage'] <- 'chloro'
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
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)

m1
###################################################################################
# formula for 2 state model with transition probabilities 
#formula <- ~ bath 
# initial parameters (obtained from nested model m1)
#Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
#m2 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
#             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
#m2
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
#formula <- ~ bath 
# initial parameters (obtained from nested model m1)
#Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
#m4 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
#             beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1)

###################################################################################
#AIC(m1,m2,m3,m4)
###################################################################################
# Chlorophyll covariate 
#stormData<-head(stormData,-5)
###################################################################################
# formula for 2 state model with transition probabilities 
#stateNames <- c("exploratory","encamped")
#formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
#Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
#m5 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
#                         beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
#m5
###################################################################################
# formula for 3 state model with transition probabilities
#stateNames <- c("transiting", "foraging", "resting")
#formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
#Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
#m6 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
#                         beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1) 
#m6
###################################################################################
#AIC(m1,m2,m3,m4,m5,m6)
###################################################################################
# Chlorophyll covariate 
#stormData<-head(stormData,-5)
###################################################################################
# formula for 2 state model with transition probabilities 
#stateNames <- c("exploratory","encamped")
#formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
#Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
#m5 <- momentuHMM::fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
#                         beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
#m5
###################################################################################
# formula for 3 state model with transition probabilities
#stateNames <- c("transiting", "foraging", "resting")
#formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
#Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
#m6 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
#                         beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1)
#m6
###################################################################################
#AIC(m1,m2,m3,m4,m5,m6,m7,m8)
###################################################################################
# decode most likely state sequence
#stormData$states <- viterbi(m3)
# derive percentage of time spent in each state
#table(stormData$states)/nrow(stormData)
#stormData$chloro <- as.numeric(stormData$chloro)


# find locations within certain lat/lon distance in r
# location of Mousa
mylat <- 60
mylon <- -1.166667
lat<-stormData$lat
lon<-stormData$lon
list1 <- data.frame(lon,lat)
list2<- data.frame(mylon,mylat)

library(geosphere)

mat <- distm(list1[,c('lon','lat')], list2[,c('mylon','mylat')], fun=distVincentyEllipsoid)
stormData$distance <- apply(mat, 1, min)/1000
head(stormData)
max(stormData$distance)

plot(stormData$x,stormData$y)
plot(stormData$x[stormData$distance>5],stormData$y[stormData$distance>5])
plot(stormData$x[stormData$distance<5],stormData$y[stormData$distance<5])

newdata <- stormData[stormData$distance >5, ]

#boxplot(newdata$chloro~newdata$states, pch = 16)
#means <- tapply(newdata$chloro,newdata$states,mean, na.rm=TRUE)
#points(means,col="red",pch=16)
#kruskal.test(newdata$chloro~newdata$states)

###################################################################################
# Fit HMMs again using the data outisde of 5km of the colony
###################################################################################
# project to UTM coordinates using package rgdal
#llcoord <- SpatialPoints(newdata[,4:5],
#                         proj4string=CRS("+proj=longlat +datum=WGS84"))
#utmcoord <- spTransform(llcoord,CRS("+proj=utm + ellps=WGS84")) # 29 = IRE or 30 = UK
## add UTM locations to data frame
#newdata$x <- attr(utmcoord,"coords")[,4]
#newdata$y <- attr(utmcoord,"coords")[,5]

#plot(newdata$x , newdata$y)

#newdata <- momentuHMM::prepData(newdata)

# remove last couple of rows from B62Blue0 which has a gigantic positive bathymetry value
newdata <- head(newdata,-2)
###################################################################################
# label states
stateNames <- c("exploratory","encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))

# fit model
m10 <- momentuHMM::fitHMM(data = newdata, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)

m10
###################################################################################
# formula for 2 state model with transition probabilities 
#formula <- ~ bath 
# initial parameters (obtained from nested model m1)
#Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
#m2 <- momentuHMM::fitHMM(data = newdata, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
#             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
#m2
###################################################################################
# label 3 states
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))

# fit model
m30 <- momentuHMM::fitHMM(data = newdata, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m30
###################################################################################
# formula for 3 state model with transition probabilities

formula <- ~ bath 
# initial parameters (obtained from nested model m1)
Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
m40 <- momentuHMM::fitHMM(data = newdata, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
                         beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1)

m40
###################################################################################
# Chlorophyll covariate
# have to take dataset out of HMM object, remove NA chloro values, then reconvert
###################################################################################
#newdata<-newdata[,c(1:11)]
newdata<-newdata[,c("lon","lat","ID", "chloro")]

# problem IDs with NAs for chlorophyll 
tail(newdata[newdata$ID=="B62Blue0",])
length(newdata$ID)
# remove NAs chlorophyll values
newdata <- newdata %>%
  dplyr:: mutate(chloro = ifelse(is.na(chloro),0,chloro))

newdata<-group_by(newdata, ID) %>%
  dplyr::mutate(first2 = min(which(chloro == 0 | row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  dplyr::select(-first2)
length(newdata$ID)
newdata<-droplevels(newdata)

newdata<-data.frame(newdata)

# make sure those NAs are gone 
tail(newdata[newdata$ID=="B62Blue0",])

# change back to momentuHMM object
# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(newdata[,1:2],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm + ellps=WGS84")) # 29 = IRE or 30 = UK
# add UTM locations to data frame
newdata$x <- attr(utmcoord,"coords")[,1]
newdata$y <- attr(utmcoord,"coords")[,2]


newdata <- momentuHMM::prepData(newdata)


# formula for 3 state model with transition probabilities
stateNames <- c("transiting", "foraging", "resting")
formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
Par0_m4 <- getPar0(model=m3, formula=formula)
# fit model
m60 <- momentuHMM::fitHMM(data = newdata, nbStates = 3, dist = dist, Par0 = Par0_m4$Par,
                         beta0=Par0_m4$beta, stateNames = stateNames, formula=formula,retryFits = 1) 
m60


########################################################################
# Compare all models 
########################################################################
AIC(m10,m30,m40, m60)