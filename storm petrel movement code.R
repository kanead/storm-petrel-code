# Storm Petrel Movement Code 

# clean everything first
rm(list=ls())

#libraries
library(adehabitatLT)
library(geosphere)
library(moveHMM)
library(rworldmap)
library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries
library(move)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
data <- read.table("allStormies.csv", header=T,sep=",")
head(data)
data<-data[,c("latitude","longitude","DateTime", "ID")]
names(data)[names(data) == 'longitude'] <- 'lon'
names(data)[names(data) == 'latitude'] <- 'lat'
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%y %H:%M", tz = "UTC")
head(data)
length(data$lat)
# remove missing data
data<-data[ ! data$lat %in% 0, ]
length(data$lat)

# plot the data
map('worldHires', c('Ireland', 'UK'),   
    xlim=c(-16,-5.5), 
    ylim=c(51,56))	
points(data$lon,data$lat,col=data$ID,pch=16, cex=0.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
       xlab="longitude",ylab="latitude")

# remove erroneous point 
data<-data[data$lat < 54.5, ]

# replot the data
map('worldHires', c('Ireland', 'UK'),   
    xlim=c(-16,-5.5), 
    ylim=c(51,56))	
points(data$lon,data$lat,col=data$ID,pch=16, cex=0.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
       xlab="longitude",ylab="latitude")

# convert into a 'move' type file 
movedata <- move(x=data$lon, y=data$lat,
             time=data$DateTime,
             data=data, proj=CRS("+proj=longlat +ellps=WGS84"), animal=data$ID)
movedata
summary(movedata)
show(movedata)
# number of relocation for each bird
n.locs(movedata)
# summary of the speed statistics in metres per second 
speedSummary(movedata)
# summary of the time statistics in hours
timeSummary(movedata, units="hours")
# summary of distance measures in metres
distanceSummary(movedata) 
# summary of angle measures in degrees
angleSummary(movedata)
# The time.lag function calculates the time lags between locations
timeLag(movedata, units="mins")

# ---------------------------------------------------------------------
# apply Hidden Markov Model to the Data 
# ---------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Interpolate the tracks so that they are measured on the same interval
# -----------------------------------------------------------------------------
# ---------------------------------------------
# try interpolating for one well behaved track 
# ---------------------------------------------
dataSample<-data[data$ID == 908, ]
head(dataSample)

# create a trajectory object using adehabitatLT
tr<-as.ltraj(data.frame(X=dataSample$lon,Y=dataSample$lat),date=dataSample$DateTime,id=dataSample$ID,typeII=T) #create trajectory
tstep<-1800 #time step we want for the interpolation, in seconds
newtr<-redisltraj(tr, u=tstep, type = "time")
head(newtr)
head(newtr[[1]])

# convert object of class ltraj to a dataframe 
df<-ld(newtr)
names(df)[names(df) == 'x'] <- 'lon'
names(df)[names(df) == 'y'] <- 'lat'
head(df)

#prepare data with moveHMM
trackData2 <- df[,c(1,2,11)]
colnames(trackData2)[3] <- c("ID")
data3 <- prepData(trackData2,type="LL",coordNames=c("lon","lat"))
plot(data3,compact=T)

#apply two state HMM
## initial parameters for gamma and von Mises distributions
mu0 <- c(1,4) # step mean (two parameters: one for each state)
sigma0 <- c(0.5,1) # step SD
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(0.7,1.5) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate

m1
plot(m1)

# -----------------------------------------------------------------------------
# Currently this does not work for all tracks being interpolated 
# -----------------------------------------------------------------------------
idx = c(900,902,908,910)
dataSample2<-data[data$ID %in% idx,] 
head(dataSample2)
tail(dataSample2)

# create a trajectory object using adehabitatLT
tr<-as.ltraj(data.frame(X=dataSample2$lon,Y=dataSample2$lat),date=dataSample2$DateTime,id=dataSample2$ID,typeII=T) #create trajectory
tstep<-1800 #time step we want for the interpolation, in seconds
newtr<-redisltraj(tr, u=tstep, type = "time")
head(newtr)
head(newtr[[1]])

# convert object of class ltraj to a dataframe 
df<-ld(newtr)
names(df)[names(df) == 'x'] <- 'lon'
names(df)[names(df) == 'y'] <- 'lat'
head(df)
tail(df)

# prepare data with moveHMM
trackData2 <- df[,c(1,2,11)]
colnames(trackData2)[3] <- c("ID")
data3 <- prepData(trackData2,type="LL",coordNames=c("lon","lat"))
plot(data3,compact=T)

#apply two state HMM
# initial parameters for gamma and von Mises distributions
mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
sigma0 <- c(0.1,1) # step SD
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate

m1
plot(m1)

