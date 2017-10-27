#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)
library(data.table)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
data <- read.table("combinedData.csv", header=T,sep=",")
head(data)
data<-data[,c("Latitude","Longitude","DateTime", "BirdID","location")]
names(data)[names(data) == 'Longitude'] <- 'lon'
names(data)[names(data) == 'Latitude'] <- 'lat'
names(data)[names(data) == 'BirdID'] <- 'ID'
names(data)[names(data) == 'DateTime'] <- 'time'

# the time stamp can be a pain - set the column in excel using dd-mm-yy hh:mm:ss
data$time<-as.POSIXct(data$time, format= "%d-%m-%y %H:%M", tz = "UTC")
head(data)
length(data$lat)
# remove missing data
data<-data[ ! data$lat %in% 0, ]
length(data$lat)
# remove NAs
data<-data[complete.cases(data),]
length(data$lat)

# split up Irish data
irishdata <- data[data$location=="ireland" , ] 
irishdata<-droplevels(irishdata)

# remove erroneous point 
irishdata<-irishdata[irishdata$lat < 54.5, ]

# split up Scottish data
scottishdata <- data[data$location=="scotland" , ]
scottishdata<-droplevels(scottishdata)

# combine data back together again
# stormData<-rbind(scottishdata,irishdata)
# or select one country's tracks 
# stormData<-irishdata[irishdata$ID==908,]
# stormData<-scottishdata
 stormData<-irishdata

# count the number of relocations for each bird 
sapply(split(stormData$lat,stormData$ID),length)
# drop the birds that have fewer than x relocations 
stormData <- stormData[!(as.numeric(stormData$ID) %in% which(table(stormData$ID)<10)),]
stormData <- droplevels(stormData)
#calculate time difference in minutes
stormData <- stormData[order(stormData$ID, stormData$time),]
stormData$tdiff <- unlist(tapply(stormData$time, INDEX = stormData$ID,
                          FUN = function(x) c(0, `units<-`(diff(x), "mins"))))


# calculate the max time gap by factor level 
#maxtdiff<- tapply(stormData$tdiff, INDEX = stormData$ID,max)
#maxtdiff
#idx <- names(maxtdiff)[maxtdiff < 120]
#stormData<-stormData[stormData$ID %in% idx,]
#stormData<-droplevels(stormData)
#maxtdiff<- tapply(stormData$tdiff, INDEX = stormData$ID,max)
#maxtdiff

# split tracks that have big time gaps by first calculating the time difference
setDT(stormData)[ , ID2 := paste0(ID, cumsum(c(0, diff(time)) > 120)), by = ID]
stormData<-data.frame(stormData)
stormData$ID2<-factor(stormData$ID2)
levels(stormData$ID2)
head(stormData)
# check out the length of the tracks that have been broken up if they exceed 2.5 hours at any point 
tapply(stormData$ID2, INDEX = stormData$ID2,length)
# drop the birds that have fewer than x relocations 
stormData <- stormData[!(as.numeric(stormData$ID2) %in% which(table(stormData$ID2)<10)),]
stormData <- droplevels(stormData)
tapply(stormData$ID2, INDEX = stormData$ID2,length)
head(stormData)
maxtdiff<- tapply(stormData$tdiff, INDEX = stormData$ID2,max)
maxtdiff
#stormData <- subset(stormData, tdiff < 120)
#stormData <- droplevels(stormData)
# look at the new max time gap by factor level 
#maxtdiff<- tapply(stormData$tdiff, INDEX = stormData$ID,max)
#maxtdiff
stormData <- stormData[,c(2,1,3,5,7)]
names(stormData)[names(stormData) == 'ID2'] <- 'ID'

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(stormData[,1:2],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=29 ellps=WGS84")) # 29 = IRE or 30 = UK
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

IrishCRW <- momentuHMM::prepData(data=crwOut)

IrishCRW$datePlusOne<-IrishCRW$time + days(1)
IrishCRW$dateMinusOne<-IrishCRW$time - days(1)

IrishCRW$datePlusTwo<-IrishCRW$time + days(2)
IrishCRW$dateMinusTwo<-IrishCRW$time - days(2)

IrishCRW$datePlusThree<-IrishCRW$time + days(3)
IrishCRW$dateMinusThree<-IrishCRW$time - days(3)

IrishCRW$datePlusFour<-IrishCRW$time + days(4)
IrishCRW$dateMinusFour<-IrishCRW$time - days(4)

IrishCRW$datePlusFive<-IrishCRW$time + days(5)
IrishCRW$dateMinusFive<-IrishCRW$time - days(5)

head(IrishCRW)

write.table(IrishCRW, file="IrishCRW.csv",sep = ",")


