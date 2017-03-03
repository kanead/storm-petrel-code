library(rWind)


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

irishdata <- data[data$location=="ireland" , ] 
irishdata<-droplevels(irishdata)

# Irish birds 
map('worldHires', c('Ireland', 'UK'),   
    xlim=c(-16,-5.5), 
    ylim=c(51,56))	
points(irishdata$lon,irishdata$lat,col=irishdata$ID,pch=16, cex=0.5, map.axes(cex.axis=0.8),title("Storm Petrels"),
       xlab="longitude",ylab="latitude")

# remove erroneous point 
irishdata<-irishdata[irishdata$lat < 54.5, ]


# Load map data for visualising Wind Data

#install.packages("rworldmap")  
library(rworldmap)   
newmap <- getMap(resolution = "high")  

# If you want a time series of wind data, you can download it by using a for-in loop:  
# First, you should create an empty list where you will store all the data  

wind_serie<- list()  

# Then, you can use a wind.dl inside a for-in loop to download and store wind data of   
# the first 5 days of August 2016 at 00:00 in Atlantic region. It could take a while...  

for (d in 1:5){  
  w<-wind.dl(2016,8,d,0,-16,-5.5,51,56)  
  wind_serie[[d]]<-w  
}  


wind_serie  

# Finally, you can use wind.mean function to calculate wind average   

wind_average<-wind.mean(wind_serie)  
wind_average<-wind.fit(wind_average)  
r_average_dir<-wind2raster(wind_average, type="dir")  
r_average_speed<-wind2raster(wind_average, type="speed")  

par(mfrow=c(1,2))  

plot(r_average_dir, main="direction average")  
lines(newmap, lwd=1)  
points(0,45)
plot(r_average_speed, main="speed average")  
lines(newmap, lwd=1)

# Do the same as above but just for one time point, midnight on 23/08/17
wind_data<-wind.dl(2016,8,27,0,-16,-5.5,51,56)  

# data defined by the two vector components U and V. You can transform these data into a nicer format
wind_data<-wind.fit(wind_data)  

r_dir <- wind2raster(wind_data, type="dir")  
r_speed <- wind2raster(wind_data, type="speed")   

# Now, you can use rworldmap package to plot countries contours with your direction and speed data!  
#par(mfrow=c(1,2))  

#plot(r_dir, main="direction")  
#lines(newmap, lwd=1)  

#plot(r_speed, main="speed")  
#lines(newmap, lwd=1)  
#points(irishdata$lon[irishdata$ID=="908"],irishdata$lat[irishdata$ID=="908"],pch=16)


# Additionally, you can use arrowDir and Arrowhead (from "shape" package) functions to plot wind direction  
# over a raster graph:  

#install.packages("shape")  
library(shape)   

#dev.off()  
alpha<- arrowDir(wind_data)  
plot(r_speed, main="wind direction (arrows) and speed (colours) - Midnight 27/08/16")  
lines(newmap, lwd=1)  
Arrowhead(wind_data$lon, wind_data$lat, angle=alpha, arr.length = 0.12, arr.type="curved")  
