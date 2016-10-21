#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())

#load required libraries
library(adehabitatLT)
library(geosphere)
library(moveHMM)
library(rworldmap)
library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries
library(move)
library(RNCEP)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
data <- read.table("allStormies.csv", header=T,sep=",")
head(data)
data<-data[,c("latitude","longitude","DateTime", "ID","bathymetry")]
names(data)[names(data) == 'longitude'] <- 'lon'
names(data)[names(data) == 'latitude'] <- 'lat'
# the time stamp can be a pain - set the column in excel using dd-mm-yy hh:mm:ss
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%y %H:%M", tz = "UTC")
head(data)
length(data$lat)
# remove missing data
data<-data[ ! data$lat %in% 0, ]
length(data$lat)

#--------------------------------------------------------------------------------
# plot the data
#--------------------------------------------------------------------------------
# specify the colours
palette(c("grey","orange","blue","red","yellow","black","cyan","pink"))

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

# alternatively, plot each of the bird tracks on a separate map
mapFunc <- function(dat) {
  map('worldHires', c('Ireland', 'UK'), xlim=c(-16,-5.5), ylim=c(51,56))    
  points(dat$lon, dat$lat,pch=16, cex=0.5, map.axes(cex.axis=0.8),
         xlab="longitude",ylab="latitude")
}

op <- par(mfrow = c(2,4),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
#birdID<-as.factor(data$ID)
sapply(split(data[1:2], data$ID), mapFunc)

# plot the tracking data with bathymetry data
par(mfrow = c(1,1))
NCEP.vis.points(wx=data$bathymetry, lats=data$lat, lons=data$lon,cols=rev(heat.colors(64)),
                  title.args=list(main="Storm Petrels with Bathymetry Data"), points.args=list(cex=1),
                    image.plot.args=list(legend.args=list(text="m AMSL",adj=0, padj=-2, cex=1.15)),
                      map.args=list(xlim=c(-16,-5.5), ylim=c(51,56)))

#--------------------------------------------------------------------------------
# convert into a 'move' type file 
#--------------------------------------------------------------------------------
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

#calculate time difference in minutes
df <- data[order(data$ID, data$DateTime),]
df$tdiff <- unlist(tapply(data$DateTime, INDEX = data$ID,
                          FUN = function(x) c(0, `units<-`(diff(x), "mins"))))
#df
which(df$tdiff > 1000)
# ---------------------------------------------------------------------
# Apply Hidden Markov Model to the Data 
# ---------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Interpolate the tracks so that they are measured on the same interval
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Currently this does not work for all tracks being interpolated 
# -----------------------------------------------------------------------------
idx = c("900","902","908","910","906","906B","909","909B")
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

states <- viterbi(m1)
states[1:25]

sp <- stateProbs(m1)
head(sp)
plotStates(m1)

# compute the pseudo-residuals
pr <- pseudoRes(m1)
# time series, qq-plots, and ACF of the pseudo-residuals
plotPR(m1)
# ---------------------------------------------
# Try interpolating for one well behaved track 
# ---------------------------------------------
#dataSample<-data[data$ID == 908, ]
#head(dataSample)

# create a trajectory object using adehabitatLT
#tr<-as.ltraj(data.frame(X=dataSample$lon,Y=dataSample$lat),date=dataSample$DateTime,id=dataSample$ID,typeII=T) #create trajectory
#tstep<-1800 #time step we want for the interpolation, in seconds
#newtr<-redisltraj(tr, u=tstep, type = "time")
#head(newtr)
#head(newtr[[1]])

# convert object of class ltraj to a dataframe 
#df<-ld(newtr)
#names(df)[names(df) == 'x'] <- 'lon'
#names(df)[names(df) == 'y'] <- 'lat'
#head(df)

#prepare data with moveHMM
#trackData2 <- df[,c(1,2,11)]
#colnames(trackData2)[3] <- c("ID")
#data3 <- prepData(trackData2,type="LL",coordNames=c("lon","lat"))
#plot(data3,compact=T)

#apply two state HMM
## initial parameters for gamma and von Mises distributions
#mu0 <- c(1,4) # step mean (two parameters: one for each state)
#sigma0 <- c(0.5,1) # step SD
#stepPar0 <- c(mu0,sigma0)
#angleMean0 <- c(pi,0) # angle mean
#kappa0 <- c(0.7,1.5) # angle concentration
#anglePar0 <- c(angleMean0,kappa0)

#m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
#             formula=~1) # no covariate

#m1
#plot(m1)

#states <- viterbi(m1)
#states[1:25]

#sp <- stateProbs(m1)
#head(sp)
#plotStates(m1)



# --------------------------------------------------------------------------
# Model with an environmental covariate 
# --------------------------------------------------------------------------
# dataSample<-data[data$ID == 908, ]
# head(dataSample)

# create a trajectory object using adehabitatLT
# tr<-as.ltraj(data.frame(X=dataSample$lon,Y=dataSample$lat),date=dataSample$DateTime,id=dataSample$ID,typeII=T) #create trajectory
# tstep<-1800 #time step we want for the interpolation, in seconds
# newtr<-redisltraj(tr, u=tstep, type = "time")
# head(newtr)
# head(newtr[[1]])

# convert object of class ltraj to a dataframe 
# df<-ld(newtr)
# names(df)[names(df) == 'x'] <- 'lon'
# names(df)[names(df) == 'y'] <- 'lat'
# head(df)

# the environmental data will need to be applied to the interpolated data at this point 
# for now we'll use non interpolated data for the best track 

#prepare data with moveHMM
# trackData2 <- dataSample[,c(1,2,4,5)]
# colnames(trackData2)[3] <- c("ID")
# data3 <- prepData(trackData2,type="LL",coordNames=c("lon","lat"))
# plot(data3,compact=T)

#apply two state HMM
## initial parameters for gamma and von Mises distributions
# mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
# sigma0 <- c(0.1,1) # step SD
# stepPar0 <- c(mu0,sigma0)
# angleMean0 <- c(pi,0) # angle mean
# kappa0 <- c(1,1) # angle concentration
# anglePar0 <- c(angleMean0,kappa0)

# m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
#             formula=~1) # no covariate
# m2 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
#             formula=~bathymetry) # covariate 'bathymetry'

## Model selection using the AIC
# print(AIC(m1,m2))

# m1
# plot(m1)
# m2
# plot(m2)

