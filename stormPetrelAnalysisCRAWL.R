#############################################################################
# STORM PETREL ANALYSIS
#############################################################################
# Clean the data
# Subset for Scotland
# Subset for median time difference between relocations < 30 mins 
# Split tracks with big time gaps 
# Remove tracks with < 20 relocations 
# Run the CRAWL model 
# Interpolate to 20 min relocations using adehabitatLT
# Export for use in Movebank 
#############################################################################
rm(list=ls())
#############################################################################
# load data
#############################################################################
setwd("C:/Users/akane/Desktop/Science/Manuscripts/Storm Petrels/Tracking Data")
mydata<-read.csv("combinedData.csv",header = T,sep = ",")
head(mydata)
#############################################################################
# drop unused columns 
#############################################################################
drops <- c("TagID","Stage","Day","Month","Year","Hour","Minute","Year","Second","Date","Timne","TripNumber","FixInterval",
           "FixIntervalHrs","ElapsedSecondsPerDay","NumberSatellites","Departure.time")
mydata<-mydata[ , !(names(mydata) %in% drops)]
head(mydata)
#############################################################################
# Convert DateTime to the proper time class
#############################################################################
mydata$DateTime <- as.POSIXct(strptime(mydata$DateTime, "%d/%m/%Y %H:%M"), "GMT") 
head(mydata)
#############################################################################
# Focus on Scottish data only 
#############################################################################
mydata <- mydata[mydata$location=="scotland",]
mydata <- droplevels(mydata)
#############################################################################
# remove rows with no values for latitude or longitude 
#############################################################################
mydata<-mydata[complete.cases(mydata$Latitude),]
mydata<-mydata[complete.cases(mydata$Longitude),]
#############################################################################
# find locations within certain lat/lon distance in r, distance given in metres 
#############################################################################
# location of Mousa
mylat <- 60
mylon <- -1.166667
lat<-mydata$Latitude
lon<-mydata$Longitude
list1 <- data.frame(lon,lat)
list2<- data.frame(mylon,mylat)

library(geosphere)
# add a column showing the distance of the points to Mousa 
mat <- distm(list1[,c('lon','lat')], list2[,c('mylon','mylat')], fun=distVincentyEllipsoid)
mydata$distance <- apply(mat, 1, min)/1000
head(mydata)
max(mydata$distance)

plot(mydata$Longitude,mydata$Latitude)
plot(mydata$Longitude[mydata$distance>5],mydata$Latitude[mydata$distance>5])
plot(mydata$Longitude[mydata$distance<5],mydata$Latitude[mydata$distance<5])
#############################################################################
# remove points that are within 5km of Mousa 
#############################################################################
length(mydata$distance)
mydata <- mydata[mydata$distance >5, ]
length(mydata$distance)
#############################################################################
# count the number of relocations per bird 
#############################################################################
sapply(split(mydata$Latitude,mydata$BirdID),length)
length(levels(mydata$BirdID))
#############################################################################
# drop the birds that have fewer than 20 relocations 
#############################################################################
mydata <- mydata[!(as.numeric(mydata$BirdID) %in% which(table(mydata$BirdID)<20)),]
mydata <- droplevels(mydata)
length(levels(mydata$BirdID))
head(mydata)
#############################################################################
# measure the time difference between points for each bird ID using dplyr 
# - Group your data by ID
# - Compute time diffs between each timestamp in your group (the 1st time diff is NA)
#############################################################################
mydata<- mydata %>%
  group_by(BirdID) %>%
  mutate(timeDiff = c(NA, difftime(tail(DateTime, -1), head(DateTime, -1), units="min")))
#############################################################################
# find the tracks that have a median time difference of 31 mins and keep only them 
#############################################################################
tapply(mydata$timeDiff, INDEX = mydata$BirdID,median,na.rm=T)
mediantdiff<- tapply(mydata$timeDiff, INDEX = mydata$BirdID,median,na.rm=T)
subsetTracks<-which(mediantdiff<31)
subsetNames<-names(subsetTracks)
mydata<-mydata[mydata$BirdID %in% subsetNames,] 
mydata<-droplevels(mydata)
levels(mydata$BirdID)
max(mydata$timeDiff,na.rm = T)
tapply(mydata$timeDiff, INDEX = mydata$BirdID,median,na.rm=T)
#############################################################################
# Split up the tracks with median 20 min relocs if there are big time gaps 
#############################################################################
library(dplyr)
mydata<- mydata %>%
  mutate(newID = paste0(BirdID, cumsum(!is.na(timeDiff) & timeDiff > 60))) %>%
  ungroup()
mydata<-data.frame(mydata)
mydata$newID<-as.factor(mydata$newID)
head(mydata)
levels(factor(mydata$newID))
length(levels(factor(mydata$newID)))
#############################################################################
# Again drop the birds that have fewer than 20 relocations for the newID
#############################################################################
length(levels(mydata$newID))
mydata <- mydata[!(as.numeric(mydata$newID) %in% which(table(mydata$newID)<20)),]
mydata <- droplevels(mydata)
length(levels(mydata$newID))
levels(factor(mydata$newID))
#############################################################################
# Fit a crawl model 
#############################################################################
library(crawl)
library(mvtnorm)  #multivariate normal distribution
library(MASS)
library(sp)
library(rgdal)
library(foreach)
library(doParallel)

coordinates(mydata) = ~Longitude+Latitude
#define projection attributes of the raw data
proj4string(mydata) = CRS("+proj=longlat")
#now define the LAEA projection required for analysis
LAEAproj = CRS('+proj=laea +lat_0=59 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
#now apply the LAREA projection
mydata = spTransform(mydata, LAEAproj)
is.projected(mydata)
head(mydata)
plot(mydata$Longitude, mydata$Latitude)

#set initial state is a list of starting values for the mean and variance-covariance for the initial state of the model. 
#When choosing the initial parameters, it is typical to have the mean centered on the first observation with zero velocity. 
#a is the starting location for the model - the first known coordinate; 
#P is the initial velocity - a 4x4 var-cov matrix. For these data, a should correspond to the location where animal was instrumented.

# First use the displayPar function to examine the given parameters
displayPar(mov.model=~1, err.model=list(x=~1), drift.model=FALSE, fixPar=c( log(5), NA, NA), data=mydata)
ids = unique(mydata@data$newID)      #define bird IDs
#registerDoParallel(cores=2)
model_fits <-
  foreach(i = 1:length(ids)) %dopar% {
    id_data = subset(mydata,newID == ids[i])
    
    init = list(a = c(sp::coordinates(id_data)[1,1],0,
                      sp::coordinates(id_data)[1,2],0),
                P = diag(c(5 ^ 2, 1, 
                           5 ^ 2, 1)))
    
    fit <- crawl::crwMLE(
      mov.model =~ 1,
      err.model=list(x=~1),
      data = id_data,
      Time.name = "DateTime",
      initial.state = init,
      coord=c("X", "Y"),
      fixPar = c(log(5), NA,NA),
      theta=c(8,-3),
      control=list(maxit=2000,trace=1, REPORT=10),
      initialSANN=list(maxit=200, trace=1, REPORT=10), need.hess=1)
    
    fit
  }

#stopImplicitCluster()
names(model_fits) <- ids
print(model_fits)
#predict regularly spaced (in time) locations. first define start and end times, with regular interval
#This function predicts the regular-timed - hourly would be 3600 - locations along the movement path using the posterior mean and variance of the track.
#registerDoParallel(cores=2)
predData <- foreach(i = 1:length(model_fits), .combine = rbind) %dopar% {
  
  model_fits[[i]]$data$DateTime <- lubridate::with_tz(
    model_fits[[i]]$data$DateTime,"GMT")
  predTimes <- seq(min(model_fits[[i]]$data$DateTime), max(model_fits[[i]]$data$DateTime), 1800)
  tmp = crawl::crwPredict(model_fits[[i]], predTime=predTimes)
}
#stopImplicitCluster()
predData$predTimes <- intToPOSIX(predData$TimeNum)
#While a projected SpatialPointsDataFrame was passed to the crwMLE() funtion, 
#the prediction object returned from crwPredict() is not a SpatialPointsDataFrame. 
#The columns mu.x and mu.y represent the predicted coordinates which we can coerce 
#into a SpatialPointsDataFrame with the coordinates() function and then specify our projection with the proj4string function.
predData_sp <- predData
coordinates(predData_sp) <- ~mu.x+mu.y
#define projection of the predicted locs from CRAWL
LAEAproj = CRS('+proj=laea +lat_0=59 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
proj4string(predData_sp) = LAEAproj
is.projected(predData_sp)
#define longlat projection required for plotting
longlatproj =CRS('+proj=longlat +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
#now apply the Longlat projection
predData_sp = spTransform(predData_sp, longlatproj)
is.projected(predData_sp)
head(predData_sp)
#reconvert to dataframe to plotting
predData<-as.data.frame(predData_sp)
head(predData)
#############################################################################
# Plot the Crawl predicted tracks using ggplot 
#############################################################################
library(ggplot2)

theme_map = function(base_size=9, base_family="")
{
  require(grid)
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.title.x=element_text(vjust=0),
          axis.title.y=element_text(angle=90, vjust=1.25),
          axis.text.y=element_text(angle=90),
          axis.ticks=element_line(colour="black", size=0.25),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text=element_text(),
          legend.title=element_text(face="bold", hjust=0),
          panel.border=element_rect(fill=NA, colour="black"),
          panel.grid.major=element_line(colour="grey92", size=0.3, linetype=1),
          panel.grid.minor=element_blank(),
          plot.title=element_text(vjust=1),
          strip.background=element_rect(fill="grey90", colour="black", size=0.3),
          strip.text=element_text()
    )
}

p1 <- ggplot(data=predData,aes(x=mu.x,y=mu.y)) + 
  geom_path(aes(colour=newID)) + xlab("easting (meters)") +
  ylab("northing (meters)") + theme_map()
p1
#############################################################################
# Make sure its interpolated to 30 mins, CRAWL doesn't make it so  
#############################################################################
head(predData)
predData<-predData[,c("mu.x","mu.y","predTimes","newID")]
names(predData)[names(predData) == 'mu.x'] <- 'longitude'
names(predData)[names(predData) == 'mu.y'] <- 'latitude'
names(predData)[names(predData) == 'predTimes'] <- 'DateTime'
head(predData)

# create a trajectory object using adehabitatLT
library(adehabitatLT)
tr<-as.ltraj(data.frame(X=predData$longitude,Y=predData$latitude),date=predData$DateTime,id=predData$newID,typeII=T) #create trajectory
tstep<-1800 #time step we want for the interpolation, in seconds, 1800 secs = 30 mins  
newtr<-redisltraj(tr, u=tstep, type = "time")
head(newtr)
head(newtr[[1]])

# convert object of class ltraj to a dataframe 
mydata<-ld(newtr)
head(mydata)
names(mydata)[names(mydata) == 'x'] <- 'lon'
names(mydata)[names(mydata) == 'y'] <- 'lat'
names(mydata)[names(mydata) == 'date'] <- 'DateTime'
mydata<-mydata[,c("lon","lat","DateTime","id")]
head(mydata)
tail(mydata)
#############################################################################
# Export the dataframe for use in Movebank 
#############################################################################
write.csv(mydata,file = "Mousa30minCrawl.csv",row.names = F)
#############################################################################
# Read in data with Movebank Covariates 
#############################################################################
setwd("C:/Users/akane/Desktop/Science/Manuscripts/Storm Petrels/Tracking Data")
mydata<-read.csv("petrelsCrawlScotland30Mins.csv",header = T, sep = ",")
head(mydata)
mydata<-mydata[,c("location.long","location.lat","timestamp","individual.local.identifier",
                  "MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A","ETOPO1.Elevation",
                  "MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A","NASA.Distance.to.Coast")]
names(mydata)[names(mydata) == 'location.long'] <- 'lon'
names(mydata)[names(mydata) == 'location.lat'] <- 'lat'
names(mydata)[names(mydata) == 'timestamp'] <- 'DateTime'
names(mydata)[names(mydata) == 'individual.local.identifier'] <- 'ID'
names(mydata)[names(mydata) == 'MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A'] <- 'eightDayChloro'
names(mydata)[names(mydata) == 'MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A'] <- 'chloroMonth'
names(mydata)[names(mydata) == 'ETOPO1.Elevation'] <- 'bath' 
names(mydata)[names(mydata) == 'NASA.Distance.to.Coast'] <- 'coast' 
head(mydata)
sum(is.na(mydata$eightDayChloro))
sum(is.na(mydata$chloroMonth))
sum(is.na(mydata$bath))
#############################################################################
# project to UTM coordinates using package rgdal
#############################################################################
llcoord <- SpatialPoints(mydata[,1:2],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84")) # 29 = IRE or 30 = UK
# add UTM locations to data frame
mydata$x <- attr(utmcoord,"coords")[,1]
mydata$y <- attr(utmcoord,"coords")[,2]
#############################################################################
# Prepare the data to be analysed using the HMM
#############################################################################
head(mydata)
plot(mydata$x , mydata$y)
mydata<-mydata[,c("ID","x","y","chloroMonth","bath","eightDayChloro","coast")]
# delete the rows with NAs for 8 day chloro if needed 
which(is.na(mydata$eightDayChloro))
mydata <- mydata[-c(438:442), ]
sum(is.na(mydata$eightDayChloro))
# delete the row with the huge value for chlorophyll
mydata<-head(mydata,-3)
hist(mydata$eightDayChloro)
mydata <- momentuHMM::prepData(mydata)
#############################################################################
# fit HMM
# - 2 state model 
#############################################################################
stateNames <- c("transiting", "foraging")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))
# fit model
m1 <- momentuHMM::fitHMM(data = mydata, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m1
# plot(m1)
###################################################################################
# 2 state model monthly chlorophyll covariate 
###################################################################################
formula <- ~ chloroMonth 
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))
# fit model
m2 <- momentuHMM::fitHMM(data = mydata, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m2
# plot(m2)
###################################################################################
# 3 state model no covariate 
###################################################################################
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))

# fit model
m3 <- momentuHMM::fitHMM(data = mydata, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m3
# plot(m3)
###################################################################################
# 3 state model monthly chlorophyll covariate 
###################################################################################
formula <- ~ chloroMonth 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))
# fit model
m4 <- momentuHMM::fitHMM(data = mydata, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m4
# plot(m4)
###################################################################################
# 2 state model bathymetry covariate 
###################################################################################
stateNames <- c("transiting", "foraging")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ bath 
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))
# fit model
m5 <- momentuHMM::fitHMM(data = mydata, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m5
# plot(m5)
###################################################################################
# 3 state model bathymetry covariate 
###################################################################################
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ bath 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))
# fit model
m6 <- momentuHMM::fitHMM(data = mydata, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m6
# plot(m6)
###################################################################################
# 2 state model 8 day chloro covariate 
###################################################################################
stateNames <- c("transiting", "foraging")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ eightDayChloro 
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))
# fit model
m7 <- momentuHMM::fitHMM(data = mydata, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m7
# plot(m7)
###################################################################################
# 3 state model 8 day chloro covariate 
###################################################################################
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ eightDayChloro 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))
# fit model
m8 <- momentuHMM::fitHMM(data = mydata, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m8
# plot(m8)

###################################################################################
# 4 state model no covariate 
###################################################################################
stateNames <- c("transiting","foraging","scoping","resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m3 <- list(step=c(20000,10000,5000,1000,1000,5000,500,500),angle=c(10,5,5,1))
# fit model
m9 <- momentuHMM::fitHMM(data = mydata, nbStates = 4, dist = dist, Par0 = Par0_m3,
                          estAngleMean = list(angle=FALSE), stateNames = stateNames,retryFits = 1)
m9
# plot(m9)
###################################################################################
# 4 state model monthly chloro covariate 
###################################################################################
stateNames <- c("transiting","foraging","scoping","resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ chloroMonth 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,5000,1000,1000,5000,500,500),angle=c(10,5,5,1))
# fit model
m10 <- momentuHMM::fitHMM(data = mydata, nbStates = 4, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m10
# plot(10)
###################################################################################
# 4 state model 8 day chloro covariate 
###################################################################################
stateNames <- c("transiting","foraging","scoping","resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ eightDayChloro 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,5000,1000,1000,5000,500,500),angle=c(10,5,5,1))
# fit model
m11 <- momentuHMM::fitHMM(data = mydata, nbStates = 4, dist = dist, Par0 = Par0_m3,
                          estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m11
# plot(m11)
###################################################################################
# 4 state model bathymetry covariate 
###################################################################################
stateNames <- c("transiting","foraging","scoping","resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ bath 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,5000,1000,1000,5000,500,500),angle=c(10,5,5,1))
# fit model
m12 <- momentuHMM::fitHMM(data = mydata, nbStates = 4, dist = dist, Par0 = Par0_m3,
                          estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m12
# plot(m12)




###################################################################################
# 3 state model 8 day chloro covariate 
###################################################################################
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# formula
formula <- ~ coast 
# initial parameters
Par0_m3 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))
# fit model
m13 <- momentuHMM::fitHMM(data = mydata, nbStates = 3, dist = dist, Par0 = Par0_m3,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames,formula = formula, retryFits = 1)
m13
# plot(m13)


###################################################################################
# Compare models using AIC
###################################################################################
AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)
