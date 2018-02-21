#######################################
# Increment time on Storm Petrel tracks
#######################################
rm(list=ls())
library(momentuHMM)
library(lubridate)
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\crawl SP-Movebank")
stormData <- read.csv("crawl SP-3493787305078051578.csv",header = T , sep = ",")
head(stormData)

stormData<-stormData[,c("location.long","location.lat","timestamp", "individual.local.identifier","ETOPO1.Elevation",
                        "NASA.Distance.to.Coast", "MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A")]
names(stormData)[names(stormData) == 'location.long'] <- 'lon'
names(stormData)[names(stormData) == 'location.lat'] <- 'lat'
names(stormData)[names(stormData) == 'individual.local.identifier'] <- 'ID'
names(stormData)[names(stormData) == 'ETOPO1.Elevation'] <- 'bath'
names(stormData)[names(stormData) == 'NASA.Distance.to.Coast'] <- 'coast'
names(stormData)[names(stormData) == 'MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A'] <- 'chloro'
stormData$timestamp<-as.POSIXct(stormData$timestamp, format= "%Y-%m-%d %H:%M", tz = "UTC")
head(stormData)

stormData$datePlusOne<-stormData$timestamp + days(1)
stormData$dateMinusOne<-stormData$timestamp - days(1)

stormData$datePlusTwo<-stormData$timestamp + days(2)
stormData$dateMinusTwo<-stormData$timestamp - days(2)

stormData$datePlusThree<-stormData$timestamp + days(3)
stormData$dateMinusThree<-stormData$timestamp - days(3)

stormData$datePlusFour<-stormData$timestamp + days(4)
stormData$dateMinusFour<-stormData$timestamp - days(4)

stormData$datePlusFive<-stormData$timestamp + days(5)
stormData$dateMinusFive<-stormData$timestamp - days(5)

head(stormData)

write.csv(stormData, file="incrementDateCrawlSP.csv")
