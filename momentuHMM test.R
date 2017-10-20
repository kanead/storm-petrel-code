#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)
library(adehabitatLT)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking data")
stormData <- read.table("908test.csv", header=T,sep=",")
head(stormData)
stormData<-stormData[,c("Latitude","Longitude","DateTime", "BirdID","location")]
names(stormData)[names(stormData) == 'Longitude'] <- 'lon'
names(stormData)[names(stormData) == 'Latitude'] <- 'lat'
names(stormData)[names(stormData) == 'BirdID'] <- 'ID'
names(stormData)[names(stormData) == 'DateTime'] <- 'time'

stormData$time<-as.POSIXct(stormData$time, format= "%d/%m/%Y %H:%M", tz = "UTC")
#head(stormData)
length(stormData$lat)
stormData<-stormData[complete.cases(stormData[,1:2 ]),]    
length(stormData$lat)

stormData <- stormData[,c(2,1,3,4)]

# create a trajectory object using adehabitatLT
tr<-as.ltraj(data.frame(X=stormData$lon,Y=stormData$lat),date=stormData$time,id=stormData$ID,typeII=T) #create trajectory
tstep<-1800 #time step we want for the interpolation, in seconds, 1800 secs = 30 mins 
newtr<-redisltraj(tr, u=tstep, type = "time")
#head(newtr)
#head(newtr[[1]])

# convert object of class ltraj to a dataframe 
df<-ld(newtr)
names(df)[names(df) == 'x'] <- 'lon'
names(df)[names(df) == 'y'] <- 'lat'
#head(df)

# prepare data with moveHMM
trackData <- df[,c(1,2,3)]
#head(trackData)
trackData <- prepData(trackData,type="LL",coordNames=c("lon","lat"))

# label states
stateNames <- c("exploratory","encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(10,5,1,2),angle=c(10,5)) # # it's mean1, mean2, sd1, sd2 for step lengths

# fit model
m1 <- fitHMM(data =  trackData, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 10)
m1
plot(m1)
 