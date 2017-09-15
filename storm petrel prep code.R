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
library(data.table) # for renaming tracks with big time gaps 
library(dplyr)
# make sure to use detach(package:plyr)


setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
data <- read.table("combinedData.csv", header=T,sep=",")
head(data)
data<-data[,c("Latitude","Longitude","DateTime", "BirdID","location")]
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
# remove NAs
data<-data[complete.cases(data),]
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
# split up Irish data
irishdata <- data[data$location=="ireland" , ] 
irishdata<-droplevels(irishdata)

# remove erroneous point 
irishdata<-irishdata[irishdata$lat < 54.5, ]

# find max distance from colony for each of the tracks for the Irish birds 
# High Island coordinates 53.5464, -10.2572	
# there are 2 other functions to measure spatial distance in the geosphere package but they give the same values
myfuncCosineIre<-function(x){max(distCosine(c(-10.2572,53.5464), cbind(x$lon, x$lat)))/1000}
sapply(split(irishdata[2:1],irishdata$ID), myfuncCosineIre)

# split up Scottish data
scottishdata <- data[data$location=="scotland" , ]
scottishdata<-droplevels(scottishdata)

# find max distance from colony for each of the tracks for the Scottish birds 
# Mousa coordinates 60, -1.166667	
myfuncCosineScot<-function(x){max(distCosine(c(-1.166667,60), cbind(x$lon, x$lat)))/1000}
sapply(split(scottishdata[2:1],scottishdata$ID), myfuncCosineScot)

# combine data back together again
data<-rbind(scottishdata,irishdata)
head(data)
tail(data)
length(data$lat)
data<-data[complete.cases(data),]
length(data$lat)

# count the number of relocations for each bird 
sapply(split(data$lat,data$ID),length)
# count the number of relocations for each bird 
median(sapply(split(data$lat,data$ID),length))
# drop the birds that have fewer than x relocations 
data <- data[!(as.numeric(data$ID) %in% which(table(data$ID)<10)),]
data <- droplevels(data)

#calculate time difference in minutes
df <- data[order(data$ID, data$DateTime),]
df$tdiff <- unlist(tapply(data$DateTime, INDEX = data$ID,
                          FUN = function(x) c(0, `units<-`(diff(x), "mins"))))
df
length(df$lat)
length(which(df$tdiff > 60))
# calculate the median temporal resolution by factor level 
mediantdiff<- tapply(df$tdiff, INDEX = df$ID,median)
# calculate the max temporal resolution by factor level 
maxtdiff<- tapply(df$tdiff, INDEX = df$ID,max)

# split tracks that have big time gaps by first calculating the time difference
setDT(data)[ , ID2 := paste0(ID, cumsum(c(0, diff(DateTime)) > 240)), by = ID]
data<-data.frame(data)
data$ID2<-factor(data$ID2)
levels(data$ID2)
head(data)
# check out the length of the tracks that have been broken up if they exceed 5 hours at any point 
tapply(data$ID2, INDEX = data$ID2,length)

# drop the newly defined tracks that have fewer than x relocations 
data <- data[!(as.numeric(data$ID2) %in% which(table(data$ID2)<10)),]
data <- droplevels(data)

tapply(data$ID2, INDEX = data$ID2,length)

data<-data[data$location=="scotland", ]
data <- data[,c(1,2,3,6)]
names(data)[names(data) == "ID2"] <- "ID"
head(data)
data <- data[,c(4,3,2,1)]
head(data)
data <- droplevels(data)

# create a trajectory object using adehabitatLT
tr<-as.ltraj(data.frame(X=data$lon,Y=data$lat),date=data$DateTime,id=data$ID,typeII=T) #create trajectory
tstep<-3600 #time step we want for the interpolation, in seconds, 3600 secs = 1 hour 
newtr<-redisltraj(tr, u=tstep, type = "time")
head(newtr)
head(newtr[[1]])

# convert object of class ltraj to a dataframe 
df<-ld(newtr)
names(df)[names(df) == 'x'] <- 'lon'
names(df)[names(df) == 'y'] <- 'lat'
head(df)
tail(df)

# can export this dataframe and use it to get remote sensing data
write.table(df, file = "subsetTracks_1Hr_Interpolation.csv", row.names=F, sep=",")




write.csv(stormDataScot,"scottishStormPetrels.csv")

# B04Blue1 has suspect convergence so it's dropped
length(stormDataScot$lat)
stormDataScot<-stormDataScot[!(stormDataScot$ID=="B04Blue1"), ]
stormDataScot <- droplevels(stormDataScot)
length(stormDataScot$lat)

# keep only B2Pink0 
stormDataScot<-stormDataScot[stormDataScot$ID=="B2Pink0", ]
stormDataScot <- droplevels(stormDataScot)

# convert times from factors to POSIX
stormDataScot$DateTime <- as.POSIXct(stormDataScot$DateTime,tz="GMT")
names(stormDataScot)[names(stormDataScot) == "DateTime"] <- "time"

# project to UTM coordinates using package rgdal
library(rgdal)
llcoord <- SpatialPoints(stormDataScot[,3:4],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))
# add UTM locations to data frame
stormDataScot$x <- attr(utmcoord,"coords")[,1]
stormDataScot$y <- attr(utmcoord,"coords")[,2]

# initial parameters for crawl fit
inits <- list(a = c(stormDataScot$x[1],0,stormDataScot$y[1],0),
              P = diag(c(5000^2, 10*3600^2, 5000^2, 10*3600^2)))
# fit crawl model
crwOut <- crawlWrap(obsData=stormDataScot, timeStep="hour", initial.state=inits,
                    theta=c(4,-10), fixPar=c(NA,NA))


testData <- prepData(data=crwOut, covNames="locType")


 