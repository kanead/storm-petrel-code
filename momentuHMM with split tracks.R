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
median(sapply(split(data$lat,data$ID),length))
# drop the birds that have fewer than x relocations 
data <- data[!(as.numeric(data$ID) %in% which(table(data$ID)<20)),]
data <- droplevels(data)

#calculate time difference in minutes
df <- data[order(data$ID, data$DateTime),]
df$tdiff <- unlist(tapply(data$DateTime, INDEX = data$ID,
                          FUN = function(x) c(0, `units<-`(diff(x), "mins"))))

# which(df$tdiff > 1000)
# calculate the median temporal resolution by factor level 
mediantdiff<- tapply(df$tdiff, INDEX = df$ID,median)
maxdiff<- tapply(df$tdiff, INDEX = df$ID,max)
df2 = ifelse(mediantdiff > 40 ,0, 1)
sort(df2)

# tracks with median temporal resolution less than 40 mins 
idx<-c("B04Blue", "B232Pink", "B46Pink", "B50Blue", "B53Pink","B62Blue","B65Blue","B72Pink","B73Blue","B79Blue", 
       "N237Pink","N62Blue","900","902","906","908","909","910")

# remove the tracks with a temporal resolution > 40 mins as set out with idx
data<-data[data$ID %in% idx,] 
data<-droplevels(data)

# split tracks that have big time gaps by first calculating the time difference
setDT(data)[ , ID2 := paste0(ID, cumsum(c(0, diff(DateTime)) > 150)), by = ID]
data<-data.frame(data)
data$ID2<-factor(data$ID2)
levels(data$ID2)
head(data)
# check out the length of the tracks that have been broken up if they exceed 2.5 hours at any point 
tapply(data$ID2, INDEX = data$ID2,length)

# remove the tracks that are still misbehaving. First, specify them. They have one reloc by itself or are split into multiple tiny tracks
idx2<-c("9021","9091","9092","9093","9094","9095","9096","B2Pink2", "B46Pink2","B46Pink3","B53Pink3","B53Pink4","B72Pink3","B73Blue4","N237Pink3")

# write a function that is the opposite of the match function %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# We remove them using the new function and drop their levels again here 
data<-data[data$ID2 %not in% idx2,] 
data<-droplevels(data)

head(data)
stormData<-data
stormData <- stormData[,c(2,1,3,5,6)]
names(stormData)[names(stormData) == 'ID2'] <- 'ID'
names(stormData)[names(stormData) == 'DateTime'] <- 'time'
sapply(split(stormData$lat,stormData$ID),length)
#stormData <- stormData[stormData$location=="ireland" , ] 
stormData <- stormData[stormData$location=="scotland" , ] 
stormData<-droplevels(stormData)

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(stormData[,1:2],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm + ellps=WGS84")) # 29 = IRE or 30 = UK
# add UTM locations to data frame
stormData$x <- attr(utmcoord,"coords")[,1]
stormData$y <- attr(utmcoord,"coords")[,2]

plot(stormData$x , stormData$y)


# drop problematic tracks that are giving NaN variance estimates with crawl  
#idx<-c("B2Pink", "B232bBlue", "B35Pink", "B232Pink", "B46Pink", "B76Blue", "B79Blue", "910")
#`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# We remove them using the new function and drop their levels again here 
#stormData<-stormData[stormData$ID %not in% idx,] 
stormData<-droplevels(stormData)

# initial parameters for crawl fit
inits <- list(a = c(stormData$x[1],0,stormData$y[1],0),
              P = diag(c(5000^2, 10*3600^2, 5000^2, 10*3600^2)))
# fit crawl model
# theta is composed of 2 parameters:
# sigma which controls the variability in velocity and 
# beta which is the autocorrelation parameter 
crwOut <- crawlWrap(obsData=stormData, timeStep="30 mins", initial.state=inits,
                    theta=c(10,-4), fixPar=c(NA,NA), retryFits = 10)

stormData <- momentuHMM::prepData(data=crwOut)

stormData$hour <- as.integer(strftime(stormData$time, format = "%H", tz="GMT"))
acf(stormData$step[!is.na(stormData$step)],lag.max=30)

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

# label 3 states
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m2 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))

# fit model
m2 <- momentuHMM::fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m2,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)

m2
