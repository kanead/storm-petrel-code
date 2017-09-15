#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)


setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\Storm Petrel Scottish Data-4040838653279469394")

data <- read.table("Storm Petrel Scottish Data-4040838653279469394.csv", header=T,sep=",")
data$timestamp<-as.POSIXct(data$timestamp, format= "%Y-%m-%d %H:%M", tz = "UTC")
head(data)

stormDataScot <- data[,c(3,4,5,8,11)]
names(stormDataScot)[names(stormDataScot) == "tag.local.identifier"] <- "ID"
names(stormDataScot)[names(stormDataScot) == "timestamp"] <- "time"
names(stormDataScot)[names(stormDataScot) == "MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A"] <- "chloro"
names(stormDataScot)[names(stormDataScot) == "location.long"] <- "lon"
names(stormDataScot)[names(stormDataScot) == "location.lat"] <- "lat"
head(stormDataScot)
stormDataScot <- stormDataScot[,c(4,1,2,3,5)]
head(stormDataScot)

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(stormDataScot[,3:4],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))
# add UTM locations to data frame
stormDataScot$x <- attr(utmcoord,"coords")[,1]
stormDataScot$y <- attr(utmcoord,"coords")[,2]

# keep only one track to test
stormDataScot<-stormDataScot[stormDataScot$ID=="N237Pink2", ]
stormDataScot <- droplevels(stormDataScot)

# initial parameters for crawl fit
inits <- list(a = c(stormDataScot$x[1],0,stormDataScot$y[1],0),
              P = diag(c(5000^2, 10*3600^2, 5000^2, 10*3600^2)))
# fit crawl model
crwOut <- crawlWrap(obsData=stormDataScot, timeStep="30 mins", initial.state=inits,
                    theta=c(4,-10), fixPar=c(NA,NA))

stormData <- prepData(data=crwOut, covNames="chloro")

stormData$hour <- as.integer(strftime(stormData$time, format = "%H", tz="GMT"))
acf(stormData$step[!is.na(stormData$step)],lag.max=30)

# label states
stateNames <- c("exploratory","encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")
# initial parameters
Par0_m1 <- list(step=c(10,5,1,2),angle=c(0.7,0.3))

# fit model
m1 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames)
