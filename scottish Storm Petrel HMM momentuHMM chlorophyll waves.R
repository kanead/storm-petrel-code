#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)
library(dplyr)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation-Waves")
df<-read.csv("subsetTracks30MinInterpolation-waves.csv",header=T,sep=",")
names(df)

# rename columns
names(df)[names(df) == 'ETOPO1.Elevation'] <- 'bath'
names(df)[names(df) == 'ECMWF.Interim.Full.Daily.SFC.Mean.Wave.Direction'] <- 'wave'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'
names(df)[names(df) == 'tag.local.identifier'] <- 'ID'
length(df$ID)

# remove positive bathymetry values
df<-group_by(df, ID) %>%
  mutate(first2 = min(which(bath > 0 | row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  select(-first2)
length(df$ID)

sapply(split(df$lat,df$ID),length)
# drop the birds that have fewer than x relocations 
df <- df[!(as.numeric(df$ID) %in% which(table(df$ID)<2)),]
df <- droplevels(df)
length(df$ID)

# split up data by location - Ireland or Scotland 
df<-data.frame(df)
df <- df[df$lon>-3 , ] 
df<-droplevels(df)

# convert times from factors to POSIX
head(df)
df$date<-as.POSIXct(df$date, format= "%d/%m/%Y %H:%M", tz = "GMT")
head(df)

# prepare data with moveHMM
stormData <- df[,c(4,5,8,13,15)]
head(stormData)
stormData <- prepData(stormData,type="LL",coordNames=c("lon","lat"))
# plot(stormData,compact=T)

# add cosinor covariate based on hour of day
stormData$hour <- as.integer(strftime(stormData$date, format = "%H", tz="GMT"))
head(stormData)

# label 2 states
stateNames <- c("exploratory", "encamped")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(10,5,1,2),angle=c(10,5)) # it's mean1,mean2,sd1,sd2 for step lengths

# fit model
m1 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)


# formula for transition probabilities without time 
formula <- ~ wave 
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)


# formula for transition probabilities with time 
formula <- ~ wave * cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m3 <- getPar0(model=m1, formula=formula)
# fit model
m3 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0=Par0_m3$beta, stateNames = stateNames, formula=formula,retryFits = 1)

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ wave * cosinor(hour, period = 24),
                       sd = ~ wave * cosinor(hour, period = 24)),
           angle = list(concentration = ~ wave))
# initial parameters (obtained from nested model m3)
Par0_m4 <- getPar0(model=m3, formula=formula, DM=DM)

# fit model
m4 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m4$Par,
             beta0 = Par0_m4$beta, DM = DM, stateNames = stateNames,
             formula = formula)

# label 3 states
stateNames <- c("exploratory", "encamped", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m5 <- list(step=c(10,5,1,1,2,1),angle=c(10,5,1))

# fit model
m5 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m5,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)



# formula for transition probabilities without time 
formula <- ~ wave 
# initial parameters (obtained from nested model m1)
Par0_m6 <- getPar0(model=m5, formula=formula)
# fit model
m6 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m6$Par,
             beta0=Par0_m6$beta, stateNames = stateNames, formula=formula,retryFits = 1)


# formula for transition probabilities with time 
formula <- ~ wave * cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m7 <- getPar0(model=m5, formula=formula)
# fit model
m7 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m7$Par,
             beta0=Par0_m7$beta, stateNames = stateNames, formula=formula,retryFits = 1)

# plot(m6, plotCI = FALSE, covs = data.frame(hour=7))

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ wave * cosinor(hour, period = 24),
                       sd = ~ wave * cosinor(hour, period = 24)),
           angle = list(concentration = ~ wave))
# initial parameters (obtained from nested model m6)
Par0_m8 <- getPar0(model=m6, formula=formula, DM=DM)

# fit model
m8 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m8$Par,
             beta0 = Par0_m8$beta, DM = DM, stateNames = stateNames,
             formula = formula)

# plot(m8, plotCI = TRUE, covs = data.frame(hour=10))
AIC(m1,m2,m3,m4,m5,m6,m7,m8)

plot(m8,covs = data.frame(hour=6))
