# MarMap Package
# http://www.molecularecologist.com/2015/07/marmap/

# clean everything first
rm(list=ls())

# Load package
library(marmap)

# load tracking data

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
# select the Irish data only 
irishdata <- data[data$location=="ireland" , ] 
irishdata<-droplevels(irishdata)
irishdata<-irishdata[irishdata$lat < 54.5, ]

# rename it 
track1<-irishdata
# match it to example data
myvars <- c("lon", "lat","ID")
track1 <- track1[myvars]
head(track1)

# Download bathymetric data and save on disk
bat <- getNOAA.bathy(-16,-5.5, 51, 56, res = 1, keep = TRUE)

# Get depth profile along both tracks and remove positive depths (since random fake values can be on land)
# path.profile() gets depth value for all cells of the bathymetric grid below the gps tracks
path1 <- path.profile(track1[,-3], bat) ; path1 <- path1[-path1[,4]>0,]

# Get depth values for each gps tracking point
# get.depth() retrieve depth values only for gps tracking points
depth1 <- get.depth(bat, track1$lon, track1$lat, locator = FALSE, distance = TRUE) 

# Add depth values to tracks 1 and 2 and remove positive depths
track1$depth <- depth1$depth ; track1 <- track1[-track1$depth > 0,]

# Plot
# layout(matrix(c(1, 2, 1, 3, 4, 4), ncol = 2, byrow = TRUE), height = c(1, 1, 1))

# Select colours
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

## Bathymetric map with gps tracks
plot(bat, land = TRUE, image = TRUE, lwd = 0.2, bpal = list(c(min(bat,na.rm=TRUE), 0, blues), c(0, max(bat, na.rm=TRUE), greys)), ylim=c(51,56))
plot(bat, deep = 0, shallow = 0, step = 0, lwd = 0.8, add = TRUE)
lines(track1[track1$ID=="900",], col = "purple",lwd = 1.5 )
lines(track1[track1$ID=="902",], col = "red" ,lwd = 1.5 )
lines(track1[track1$ID=="906",], col = "blue",lwd = 1.5 )
lines(track1[track1$ID=="908",], col = "green4",lwd = 1.5 )
lines(track1[track1$ID=="909",], col = "orange",lwd = 1.5 )
lines(track1[track1$ID=="910",], col = "cyan",lwd = 1.5 )
legend("topright", legend = c("900", "902", "906", "908", "909", "910"), lwd = 1, col = c("purple", "red", "blue", "green4", "orange", "yellow4"), pch = 1, pt.cex = 0.5, bg="white")

## Depths profiles along gps tracks
plotProfile(path1, main = "Track 1")

## Depth as a function of time since the deployment of gps tags
par(mar=c(5,4,1,2))
plot(track1$time, track1$depth, xlab = "Time (arbitrary unit)", ylab = "Depth (m)", type = "o", cex = 0.5, col = "brown3", xlim = range(c(track1$time, track2$time)), ylim = range(c(track1$depth, track2$depth)))
lines(track2$time, track2$depth, type = "o", cex = 0.5, col = "darkorange1")
legend("bottomright", legend = c("Track 1", "Track 2"), lwd = 1, col = c("brown3", "darkorange1"), pch = 1, pt.cex = 0.5)