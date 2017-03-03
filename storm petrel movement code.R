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
sapply(split(irishdata[2:1],irishdata$ID), myfuncCosine)

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
df
# which(df$tdiff > 1000)
# calculate the median temporal resolution by factor level 
mediantdiff<- tapply(df$tdiff, INDEX = df$ID,median)

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
idx2<-c("9021","9091","9092","9093","9094","9095","9096","B46Pink2","B46Pink3","B53Pink3","B53Pink4","B72Pink3","B73Blue4","N237Pink3")

# write a function that is the opposite of the match function %in%
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# We remove them using the new function and drop their levels again here 
data<-data[data$ID2 %not in% idx2,] 
data<-droplevels(data)

# create a trajectory object using adehabitatLT
tr<-as.ltraj(data.frame(X=data$lon,Y=data$lat),date=data$DateTime,id=data$ID2,typeII=T) #create trajectory
tstep<-1800 #time step we want for the interpolation, in seconds, 1800 secs = 30 mins 
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
write.table(df, file = "subsetTracks30MinInterpolation.csv", row.names=F, sep=",")

# ----------------------------------------------------------------------------------------- 
# read dataframe back in with remote sensing data appended from Movebank for chlorophyll
# ----------------------------------------------------------------------------------------
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation- 8 Day Chlorophyll")
df<-read.csv("subsetTracks30MinInterpolation-8770236889496130104.csv",header=T,sep=",")

# rename columns
names(df)[names(df) == 'MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A'] <- 'chloro'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'
names(df)[names(df) == 'tag.local.identifier'] <- 'ID'

# prepare data with moveHMM
trackData2 <- df[,c(4,5,8,11)]
head(trackData2)
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

# fit a 2 state model with no covariate 
m1 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~1) 

# fit a 2 state model with covariate included
m2 <- fitHMM(data=data3,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~chloro) 

# look at the output of each model and plot it
m1
plot(m1)

m2
plot(m2)

## Model selection using the AIC
print(AIC(m1,m2))

#initial parameters for distributions (model with 3 states)
mu0 <- c(0.1,1,1)                     #step mean
sigma0 <- c(0.1,1,1)                  #step SD
#zeromass0 <- c(0.2,0.1,0.1)           #step zero-mass
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0,pi)              #angle mean
kappa0 <- c(1,2,2)                    #angle concentration
anglePar0 <- c(angleMean0,kappa0)

#HMM fitting 3 states
m3 <- fitHMM(data=data3,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0)

# fit a 3 state model with covariate included
m4 <- fitHMM(data=data3,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~chloro) 

# look at the output of each model and plot it
m3
plot(m3)
m4
plot(m4)

# compare AIC scores of the 4 models 
print(AIC(m1,m2,m3,m4))

# save the viterbi states to the interpolated track 
statesm4 <- viterbi(m4)
data3$viterbi <- statesm4
head(data3)

# export the interpolated data with the Viterbi sequence 
trackPlusViterbi <- data.frame(data3)
write.table(trackPlusViterbi, file = "tracksPlusViterbi.csv", row.names=F, sep=",")

# load this data back in here if necessary 
# data3 <- read.table("tracksPlusViterbi.csv", header=T,sep=",")

# plot the track with the viterbi sequence 
plot(data3$x[data3$ID=="9080"],data3$y[data3$ID=="9080"], col=c("goldenrod3", "lightskyblue", "seagreen3")[data3$viterbi[data3$ID=="9080"]],
     pch=16,lty=1,xlab="longitude", ylab="latitude")

# plot the track with the viterbi sequence on a map 
# create the map 
library(rworldmap)   
newmap <- getMap(resolution = "high")  
plot(newmap,
       xlim = c(-15.5, -8),
       ylim = c(51.5, 55.5)
     )

# add the points of the track 
points(data3$x[data3$ID=="9080"],data3$y[data3$ID=="9080"], col=c("goldenrod3", "lightskyblue", "seagreen3")[data3$viterbi[data3$ID=="9080"]],
     pch=20,lty=1,xlab="longitude", ylab="latitude")

# create a legend 
legend("topleft", # places a legend at the appropriate place 
       c("State 1","State 2", "State 3"), # puts text in the legend
       pch=20, # gives the legend appropriate symbols (lines)
       col=c("goldenrod3", "lightskyblue", "seagreen3"),
       bty='n') # gives the legend lines the correct color and width

# -------------------------------------------------------------------------------------------------
# compute the 2D binned kernel density estimate for foraging behaviour
# -------------------------------------------------------------------------------------------------
library(adehabitatHR)
data3State3 <- data3[data3$viterbi=="3" , ] 
data3State3<-droplevels(data3State3)

spdf <- SpatialPointsDataFrame(coordinates(cbind(data3State3$x, data3State3$y)),
                               data = data3State3)
# another quick plot
plot(spdf, pch = 19, cex = .5)


kd <- kernelUD(spdf)
image(kd)

ud <- getverticeshr(kd, percent = 90)
class(ud)
plot(ud)
# -------------------------------------------------------------------------------------------------
# Create a heat map of the behaviours of the birds over time 
# -------------------------------------------------------------------------------------------------
# here I've added the timestamp from df to data3 so I can associate the viterbi sequence with time
data3$timeStamp<-df$timestamp
head(data3)
data3$timeStamp<-as.POSIXct(data3$timeStamp)
head(data3)
# I only want the hour and minute portion of the timestamp
data3$time<-format(data3$timeStamp, "%H")
head(data3)

library("ggplot2")
library("plyr")
library("reshape2")
library("scales")

dfStates<-data3[,c("viterbi","time","ID")]

# heat map for all 3 states at the same time 
dfStates %>% 
  as.data.frame() %>% 
  group_by(viterbi, time) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(time, y = factor(viterbi), fill = n)) + 
  geom_raster() +
  scale_fill_continuous('Count', low = 'white', high = 'darkgreen', guide = 'legend') +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,seq(0, 25, 1), 24)) +
  ylab('State') + 
  xlab('Hour') +
  coord_fixed() +
  theme(legend.position = 'top',
        legend.key.height = unit(.2, 'cm'))


# separate out the states and plot them side by side 
# state 1
dfStatesSub1 <- dfStates[dfStates$viterbi=="1" , ] 
dfStatesSub1<-droplevels(dfStatesSub1)

dfStatesSub1 %>% 
  as.data.frame() %>% 
  group_by(viterbi, time) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(time, y = factor(viterbi), fill = n)) + 
  geom_raster() +
  scale_fill_continuous('Count', low = 'white', high = 'goldenrod3', guide = 'legend') +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,seq(0, 25, 1), 24)) +
  ylab('State') + 
  xlab('Hour') +
  coord_fixed() +
  theme(legend.position = 'top',
        legend.key.height = unit(.2, 'cm'))


# state 2 lightskyblue
dfStatesSub2 <- dfStates[dfStates$viterbi=="2" , ] 
dfStatesSub2<-droplevels(dfStatesSub2)

dfStatesSub2 %>% 
  as.data.frame() %>% 
  group_by(viterbi, time) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(time, y = factor(viterbi), fill = n)) + 
  geom_raster() +
  scale_fill_continuous('Count', low = 'white', high = 'lightskyblue', guide = 'legend') +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,seq(0, 25, 1), 24)) +
  ylab('State') + 
  xlab('Hour') +
  coord_fixed() +
  theme(legend.position = 'top',
        legend.key.height = unit(.2, 'cm'))


# state 3 seagreen3
dfStatesSub3 <- dfStates[dfStates$viterbi=="3" , ] 
dfStatesSub3<-droplevels(dfStatesSub3)

dfStatesSub3 %>% 
  as.data.frame() %>% 
  group_by(viterbi, time) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(time, y = factor(viterbi), fill = n)) + 
  geom_raster() +
  scale_fill_continuous('Count', low = 'white', high = 'seagreen3', guide = 'legend') +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,seq(0, 25, 1), 24)) +
  ylab('State') + 
  xlab('Hour') +
  coord_fixed() +
  theme(legend.position = 'top',
        legend.key.height = unit(.2, 'cm')) 



# ----------------------------------------------------------------------------------------- 
# read dataframe back in with remote sensing data appended from Movebank for Charnock
# ----------------------------------------------------------------------------------------
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation- Charnock")
df<-read.csv("subsetTracks30MinInterpolation-491691643874147945.csv",header=T,sep=",")

# rename columns
names(df)[names(df) == 'ECMWF.Interim.Full.Daily.SFC.FC.Charnock.Parameter'] <- 'charnock'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'
names(df)[names(df) == 'tag.local.identifier'] <- 'ID'

# prepare data with moveHMM
trackData2 <- df[,c(4,5,8,11)]
head(trackData2)
data3 <- prepData(trackData2,type="LL",coordNames=c("lon","lat"))
plot(data3,compact=T)

#initial parameters for distributions (model with 3 states)
mu0 <- c(0.1,1,1)                     #step mean
sigma0 <- c(0.1,1,1)                  #step SD
#zeromass0 <- c(0.2,0.1,0.1)           #step zero-mass
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0,pi)              #angle mean
kappa0 <- c(1,2,2)                    #angle concentration
anglePar0 <- c(angleMean0,kappa0)

#HMM fitting 3 states
m5 <- fitHMM(data=data3,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0)

# fit a 3 state model with covariate included
m6 <- fitHMM(data=data3,nbStates=3,stepPar0=stepPar0,anglePar0=anglePar0,
             formula=~charnock) 

# look at the output of each model and plot it
m5
plot(m5)
m6
plot(m6)

# compare AIC scores of the 2 models 
print(AIC(m5,m6))
 

# ----------------------------------------------------------------------------------------
#initial parameters for distributions (model with 4 states)
# ----------------------------------------------------------------------------------------
mu0 <- c(0.1,1,2,3)          #step mean
sigma0 <- c(0.01,0.1,10,1)      #step SD
#zeromass0 <- c(0.2,0.1,0.1,0.1) #step zero-mass
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(0,pi,0,pi)      #angle mean
kappa0 <- c(1,1,2,2)            #angle concentration
anglePar0 <- c(angleMean0,kappa0)

#HMM fitting 4 states
m7 <- fitHMM(data=data3,nbStates=4,stepPar0=stepPar0,anglePar0=anglePar0)
m7
plot(m7)

# compare AIC scores of the 2 models 
print(AIC(m5,m7))

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
 
 #--------------------------------------------------------------------------------
 # Plot Data 
 #--------------------------------------------------------------------------------
 # plot the Irish data
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
 
 