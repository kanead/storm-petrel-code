library(crawl)
library(mvtnorm)  #multivariate normal distribution
library(MASS)
library(sp)
library(rgdal)
library(foreach)
library(doParallel)

#install.packages("lubridate")

rm(list=ls())

#my.proj.geog = CRS('+proj=longlat +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
#my.proj.proj = CRS('+proj=laea +lat_0=-7.947 +lon_0=-14.30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data") 

Alltracks = read.csv("combinedData.csv", header=T, sep = ",")

head(Alltracks)

Alltracks$DateTime <- as.POSIXct(strptime(Alltracks$DateTime, "%d/%m/%Y %H:%M"), "GMT") 

head(Alltracks)

# Focus on Scottish data only 
Alltracks <- Alltracks[Alltracks$location=="scotland",]
Alltracks <- droplevels(Alltracks)

head(Alltracks)

# get rid of the lat and long rows with NAs
Alltracks<-Alltracks[complete.cases(Alltracks$Latitude),]
Alltracks<-Alltracks[complete.cases(Alltracks$Longitude),]


#First, tell R which columns represent the coordinates for the data

coordinates(Alltracks) = ~Longitude+Latitude

#define projection attributes of the raw data
proj4string(Alltracks) = CRS("+proj=longlat")
#now define the LAEA projection required for analysis
LAEAproj = CRS('+proj=laea +lat_0=59 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
#now apply the LAREA projection
Alltracks = spTransform(Alltracks, LAEAproj)


is.projected(Alltracks)

head(Alltracks)
plot(Alltracks$Longitude, Alltracks$Latitude)

#set initial state is a list of starting values for the mean and variance-covariance for the initial state of the model. 
#When choosing the initial parameters, it is typical to have the mean centered on the first observation with zero velocity. 
#a is the starting location for the model - the first known coordinate; 
#P is the initial velocity - a 4x4 var-cov matrix. For these data, a should correspond to the location where animal was instrumented.

# First use the displayPar function to examine the given parameters

displayPar(mov.model=~1, err.model=list(x=~1), drift.model=FALSE, fixPar=c( log(5), NA, NA), data=Alltracks)

ids = unique(Alltracks@data$BirdID)      #define bird IDs

#registerDoParallel(cores=2)

model_fits <-
  foreach(i = 1:length(ids)) %dopar% {
    id_data = subset(Alltracks,BirdID == ids[i])
    
    init = list(a = c(sp::coordinates(id_data)[1,1],0,
                      sp::coordinates(id_data)[1,2],0),
                P = diag(c(5 ^ 2, 1, 
                           5 ^ 2, 1)))
        
    fit <- crawl::crwMLE(
      mov.model =~ 1,
      err.model=list(x=~1),
      data = id_data,
      Time.name = "DateTime",
      initial.state = init,
      coord=c("X", "Y"),
      fixPar = c(log(5), NA,NA),
      theta=c(8,-3),
      control=list(maxit=2000,trace=1, REPORT=10),
      initialSANN=list(maxit=200, trace=1, REPORT=10), need.hess=1)
    
    fit
  }

#stopImplicitCluster()

names(model_fits) <- ids

print(model_fits)

#predict regularly spaced (in time) locations. first define start and end times, with regular interval

#This function predicts the regular-timed - in this case, hourly - locations along the movement path using the posterior mean and variance of the track.

#registerDoParallel(cores=2)
predData <- foreach(i = 1:length(model_fits), .combine = rbind) %dopar% {
  
  model_fits[[i]]$data$DateTime <- lubridate::with_tz(
    model_fits[[i]]$data$DateTime,"GMT")
  predTimes <- seq(min(model_fits[[i]]$data$DateTime), max(model_fits[[i]]$data$DateTime), 3600)
  tmp = crawl::crwPredict(model_fits[[i]], predTime=predTimes)
}
#stopImplicitCluster()

predData$predTimes <- intToPOSIX(predData$TimeNum)

#While a projected SpatialPointsDataFrame was passed to the crwMLE() funtion, 
#the prediction object returned from crwPredict() is not a SpatialPointsDataFrame. 
#The columns mu.x and mu.y represent the predicted coordinates which we can coerce 
#into a SpatialPointsDataFrame with the coordinates() function and then specify our projection with the proj4string function.

predData_sp <- predData
coordinates(predData_sp) <- ~mu.x+mu.y

#define projection of the predicted locs from CRAWL
LAEAproj = CRS('+proj=laea +lat_0=59 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
proj4string(predData_sp) = LAEAproj
is.projected(predData_sp)

#define longlat projection required for plotting
longlatproj =CRS('+proj=longlat +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84  +units=m +no_defs')
#now apply the Longlat projection
predData_sp = spTransform(predData_sp, longlatproj)

is.projected(predData_sp)
head(predData_sp)

#reconvert to dataframe to plotting
predData<-as.data.frame(predData_sp)

head(predData)

library(ggplot2)

theme_map = function(base_size=9, base_family="")
{
  require(grid)
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.title.x=element_text(vjust=0),
          axis.title.y=element_text(angle=90, vjust=1.25),
          axis.text.y=element_text(angle=90),
          axis.ticks=element_line(colour="black", size=0.25),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text=element_text(),
          legend.title=element_text(face="bold", hjust=0),
          panel.border=element_rect(fill=NA, colour="black"),
          panel.grid.major=element_line(colour="grey92", size=0.3, linetype=1),
          panel.grid.minor=element_blank(),
          plot.title=element_text(vjust=1),
          strip.background=element_rect(fill="grey90", colour="black", size=0.3),
          strip.text=element_text()
    )
}

p1 <- ggplot(data=predData,aes(x=mu.x,y=mu.y)) + 
  geom_path(aes(colour=BirdID)) + xlab("easting (meters)") +
  ylab("northing (meters)") + theme_map()
p1


