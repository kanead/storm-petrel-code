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
library(circular)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
data <- read.table("combinedData.csv", header=T,sep=",")
head(data)
data<-data[,c("Latitude","Longitude","DateTime", "BirdID","bathymetry","location")]
names(data)[names(data) == 'Longitude'] <- 'lon'
names(data)[names(data) == 'Latitude'] <- 'lat'
names(data)[names(data) == 'BirdID'] <- 'ID'
# the time stamp can be a pain - set the column in excel using dd-mm-yy hh:mm:ss
data$DateTime<-as.POSIXct(data$DateTime, format= "%d-%m-%y %H:%M", tz = "UTC")
head(data)
length(data$lat)
# remove missing data
data<-data[ ! data$lat %in% 0, ]
length(data$lat)


#data$lat[data$location=="ireland"] < 54.5
#y <- which(data$lat[data$location=="ireland"] < 54.5)
#head(y)
#z<-subset(data, data$location=="ireland" & data$lat < 54.5)
#length(z$lat)
#head(z)
#data<-data[data$lat[data$location=="ireland"] < 54.5, ]
#length(data$lat)
#--------------------------------------------------------------------------------
# plot the data
#--------------------------------------------------------------------------------
# specify the colours
# palette(c("grey","orange","blue","red","yellow","black","cyan","pink"))
irishdata <- data[data$location=="ireland" , ] 
irishdata<-droplevels(irishdata)

scottishdata <- data[data$location=="scotland" , ]
scottishdata<-droplevels(scottishdata)

# Irish birds 
map('worldHires', c('Ireland', 'UK'),   
    xlim=c(-16,-5.5), 
    ylim=c(51,56))	
points(irishdata$lon,irishdata$lat,col=irishdata$ID,pch=16, cex=0.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
       xlab="longitude",ylab="latitude")

# remove erroneous point 
irishdata<-irishdata[irishdata$lat < 54.5, ]

# replot the Irish data
map('worldHires', c('Ireland', 'UK'),   
    xlim=c(-16,-5.5), 
    ylim=c(51,56))	
points(irishdata$lon,irishdata$lat,col=irishdata$ID,pch=16, cex=0.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
       xlab="longitude",ylab="latitude")

# Scottish birds
map('worldHires', c('UK'),   
    xlim=c(-8,6), 
    ylim=c(56,62))	
points(scottishdata$lon,scottishdata$lat,col=scottishdata$ID,pch=16, cex=0.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
       xlab="longitude",ylab="latitude")

# stick with Irish data only for the time being 
# data <- irishdata
# data<-droplevels(data)

# alternatively, plot each of the bird tracks on a separate map

# method 1

coplot(lat ~ lon | ID, data = irishdata,pch=16)

# or 

# method 2

op <- par(mfrow = c(2,3),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
d_ply(irishdata, "ID", transform, plot(lat~lon, main = unique(ID), type = "o",pch=16))

# or 

# method 3

mapFunc <- function(data) {
  map('worldHires', c('Ireland', 'UK'), xlim=c(-16,-5.5), ylim=c(51,56))    
  points(data$lon,data$lat,pch=16, cex=.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
         xlab="longitude",ylab="latitude")
}

op <- par(mfrow = c(2,4),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

sapply(split(irishdata[2:1],irishdata$ID),mapFunc)

# or

# method 4

library(ggmap)
library(RColorBrewer)
df<-irishdata
#locate the center of the map
center<-c(mean(range(df$lon)), mean(range(df$lat)))
#in this case zoom is set by trial and error
mymap<-qmap(location = center, zoom = 2, maptype= "terrain")
mymap<-mymap + geom_point(aes(x=lon, y=lat, color=ID), data=df)
mymap<-mymap + scale_size(range = c(2, 4)) + scale_color_brewer(palette = "Set1")
mymap<-mymap + geom_path(aes(x=lon, y=lat), data=df) 

mymap<-mymap + facet_wrap(~ID, nrow =2)
print(mymap)


# plot the tracking data with bathymetry data
#par(mfrow = c(1,1))
#NCEP.vis.points(wx=data$bathymetry, lats=data$lat, lons=data$lon,cols=rev(heat.colors(64)),
#                  title.args=list(main="Storm Petrels with Bathymetry Data"), points.args=list(cex=1),
#                    image.plot.args=list(legend.args=list(text="m AMSL",adj=0, padj=-2, cex=1.15)),
#                      map.args=list(xlim=c(-16,-5.5), ylim=c(51,56)))


# combine data back together again
data<-rbind(scottishdata,irishdata)
head(data)
tail(data)
length(data$lat)
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
# idx = c("900","902","908","910","906","906B","909","909B")
#idx = c("900","902","908","910")
idx = "908"

dataSample2<-data[data$ID %in% idx,] 
dataSample2<-droplevels(dataSample2)
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

# drop the lat and long NAs
df<-df[!with(df,is.na(lat)| is.na(lon)),]

# can export this dataframe and use it to get remote sensing data
write.table(df, file = "DataInterp.csv", row.names=F, sep=",")

# read dataframe back in with remote sensing data appended 
df<-read.csv("Storm Petrel 30 Minute Interpolation-7877473479535615074.csv",header=T,sep=",")

# select one bird to test
idx = "908"
df<-df[df$tag.local.identifier %in% idx,] 

# remove last few rows where there are NAs for Chlorophyll
df<-head(df,-6)

# rename columns
names(df)[names(df) == 'MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A'] <- 'chloro'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'

# prepare data with moveHMM
trackData2 <- df[,c(4,5,11)]
head(trackData2)
#colnames(trackData2)[3] <- c("ID")
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

# alternative parameters
mu0 <- c(1,4) # step mean (two parameters: one for each state)
sigma0 <- c(0.5,1) # step SD
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(0.7,1.5) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate

m2 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~chloro) # covariate 'chlorophyll'

## Model selection using the AIC
m1
plot(m1)

m2
plot(m2)


print(AIC(m1,m2))


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
mu0 <- c(1,4) # step mean (two parameters: one for each state)
sigma0 <- c(0.5,1) # step SD
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(0.7,1.5) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

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
 mu0 <- c(0.1,1) # step mean (two parameters: one for each state)
 sigma0 <- c(0.1,1) # step SD
 stepPar0 <- c(mu0,sigma0)
 angleMean0 <- c(pi,0) # angle mean
 kappa0 <- c(1,1) # angle concentration
 anglePar0 <- c(angleMean0,kappa0)

 m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) # no covariate
 m2 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~chloro) # covariate 'chlorophyll'

## Model selection using the AIC
# print(AIC(m1,m2))

 m1
 plot(m1)
 m2
 plot(m2)

