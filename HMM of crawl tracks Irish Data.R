###################################################################################
# Storm Petrel Movement Code - Irish Data
###################################################################################
# clean everything first
rm(list=ls())
library(momentuHMM)
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\crawl SP-Movebank")
stormData <- read.csv("crawl SP Irish.csv",header = T , sep = ",")
head(stormData)

stormData<-stormData[,c("location.long","location.lat","timestamp", "individual.local.identifier","MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A", "ETOPO1.Elevation")]
names(stormData)[names(stormData) == 'location.long'] <- 'lon'
names(stormData)[names(stormData) == 'location.lat'] <- 'lat'
names(stormData)[names(stormData) == 'individual.local.identifier'] <- 'ID'
names(stormData)[names(stormData) == 'timestamp'] <- 'time'
names(stormData)[names(stormData) == 'ETOPO1.Elevation'] <- 'bath'
names(stormData)[names(stormData) == 'MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A'] <- 'chloro'

head(stormData)
stormData$time<-as.POSIXct(stormData$time, format= "%Y-%m-%d %H:%M", tz = "GMT")
stormData$hour <- as.integer(strftime(stormData$time, format = "%H", tz="GMT"))
head(stormData)

# problem IDs with NAs for chlorophyll 
tail(stormData[stormData$ID=="9020",])

# remove NAs chlorophyll values
stormData <- stormData %>%
  dplyr:: mutate(chloro = ifelse(is.na(chloro),0,chloro))

newdata<-group_by(stormData, ID) %>%
  dplyr::mutate(first2 = min(which(chloro == 0 | row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  dplyr::select(-first2)
length(stormData$ID)
stormData<-droplevels(stormData)
stormData<-stormData[stormData$chloro!=0 , ]
stormData<-data.frame(stormData)

# make sure those NAs are gone 
tail(stormData[stormData$ID=="9020",])

# project to UTM coordinates using package rgdal
llcoord <- SpatialPoints(stormData[,1:2],
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm + ellps=WGS84")) # 29 = IRE or 30 = UK
# add UTM locations to data frame
stormData$x <- attr(utmcoord,"coords")[,1]
stormData$y <- attr(utmcoord,"coords")[,2]

plot(stormData$x , stormData$y)

stormData <- momentuHMM::prepData(stormData)

# find locations within certain lat/lon distance in r, distance given in metres 
# location of High Island   
mylat <- 53.5464
mylon <- -10.2572
lat<-stormData$lat
lon<-stormData$lon
list1 <- data.frame(lon,lat)
list2<- data.frame(mylon,mylat)

library(geosphere)

mat <- distm(list1[,c('lon','lat')], list2[,c('mylon','mylat')], fun=distVincentyEllipsoid)
stormData$distance <- apply(mat, 1, min)/1000
head(stormData)
max(stormData$distance)

plot(stormData$x,stormData$y)
plot(stormData$x[stormData$distance>5],stormData$y[stormData$distance>5])
plot(stormData$x[stormData$distance<5],stormData$y[stormData$distance<5])

newdata <- stormData[stormData$distance >5, ]

###################################################################################
# Build HMM Models 
###################################################################################
# 2 state model no covariate 
###################################################################################
# label states
stateNames <- c("exploratory","encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(15000,5000,1000,5000),angle=c(10,5))

# fit model
m1 <- momentuHMM::fitHMM(data = newdata, nbStates = 2, dist = dist, Par0 = Par0_m1,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)

m1
###################################################################################
# 2 state model bathymetry covariate 
###################################################################################
formula <- ~ bath 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- momentuHMM::fitHMM(data = newdata, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
                         beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m2
###################################################################################
# 2 state model chlorophyll covariate 
###################################################################################
formula <- ~ chloro 
# initial parameters (obtained from nested model m1)
Par0_m3 <- getPar0(model=m1, formula=formula)
# fit model
m3 <- momentuHMM::fitHMM(data = newdata, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
                         beta0=Par0_m3$beta, stateNames = stateNames, formula=formula,retryFits = 1)
m3
###################################################################################
# 3 state model no covariate 
###################################################################################
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m4 <- list(step=c(20000,10000,100,1000,5000,500),angle=c(10,5,1))

# fit model
m4 <- momentuHMM::fitHMM(data = newdata, nbStates = 3, dist = dist, Par0 = Par0_m4,
                         estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)
m4
###################################################################################
# 3 state model bathymetry covariate 
###################################################################################
formula <- ~ bath 
# initial parameters (obtained from nested model m4)
Par0_m5 <- getPar0(model=m4, formula=formula)
# fit model
m5 <- momentuHMM::fitHMM(data = newdata, nbStates = 3, dist = dist, Par0 = Par0_m5$Par,
                         beta0=Par0_m5$beta, stateNames = stateNames, formula=formula,retryFits = 1)

m5
###################################################################################
# 3 state model chlorophyll covariate 
###################################################################################
formula <- ~ chloro 
# initial parameters (obtained from nested model m4)
Par0_m6 <- getPar0(model=m4, formula=formula)
# fit model
m6 <- momentuHMM::fitHMM(data = newdata, nbStates = 3, dist = dist, Par0 = Par0_m6$Par,
                         beta0=Par0_m6$beta, stateNames = stateNames, formula=formula,retryFits = 1) 
m6

########################################################################
# Compare all models 
########################################################################
AIC(m1,m2,m3,m4,m5,m6)
# decode most likely state sequence for best model
states <- viterbi(m6)
# derive percentage of time spent in each state
table(states)/nrow(newdata)
# add states to dataframe 
newdata$states<-states 

# compute pseudo-residuals for the steps and the angles
pr <- pseudoRes(m6)
# plot the ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 30)

########################################################################
# Plot the states by time 
########################################################################
library(ggplot2)
# Basic violin plot
newdata$states <- as.factor(newdata$states)
p <- ggplot(newdata, aes(x=states, y=hour, fill=states)) + 
  geom_violin(trim=TRUE, adjust = 1) + geom_boxplot(width=0.1, fill="white") +
  scale_x_discrete(labels = c('exploratory','foraging','resting'))
p + theme_set(theme_gray(base_size = 20))+
  coord_flip() + scale_fill_manual(values=c("#DAA520", "#87CEFA", "#006400"))+
  theme(legend.position="none") # Remove legend

#goldenrod3
#lightskyblue
#darkgreen
########################################################################
# Calculate home ranges for whole track and state 2 
########################################################################
library(adehabitatHR)
sp_df<-data.frame(newdata)
sp_df <- subset(sp_df, select = c(x, y))
coordinates(sp_df) <- ~x+y
plot(sp_df, axes = T)
kud <- kernelUD(sp_df, grid = 200, same4all=TRUE)
hr <- kernel.area(kud, percent = 95)
hr
image(kud)

# subset the data for state 2 alone 
sp_df2 <- newdata[newdata$states ==2, ]
sp_df2<-data.frame(sp_df2)
sp_df2 <- subset(sp_df2, select = c(x, y))
coordinates(sp_df2) <- ~x+y
plot(sp_df2, axes = T)
kud2 <- kernelUD(sp_df2, grid = 200, same4all=TRUE)
hr2 <- kernel.area(kud2, percent = 95)
hr2
image(kud2)

