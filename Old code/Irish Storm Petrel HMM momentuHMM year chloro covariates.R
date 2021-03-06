#--------------------------------------------------------------------------------
# Storm Petrel Movement Code Irish Data Monthly Chlorophyll
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)
library(dplyr)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\30minBathChloro")
df<-read.csv("subsetTracks30MinInterpolation-monthlyChloroBath.csv",header=T,sep=",")

# rename columns
names(df)[names(df) == 'ETOPO1.Elevation'] <- 'bath'
names(df)[names(df) == 'MODIS.Ocean.Aqua.OceanColor.4km.8d.Chlorophyll.A'] <- 'chloro'
names(df)[names(df) == 'MODIS.Ocean.Aqua.OceanColor.4km.Monthly.Chlorophyll.A'] <- 'chloroMonth'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'
names(df)[names(df) == 'tag.local.identifier'] <- 'ID'
length(df$ID)

# remove positive bathymetry values
df<-group_by(df, ID) %>%
  dplyr:: mutate(first2 = min(which(bath > 0 | row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  dplyr::select(-first2)
length(df$ID)


# remove NAs chlorophyll values
df <- df %>%
  dplyr:: mutate(chloro = ifelse(is.na(chloro),0,chloro))

df<-group_by(df, ID) %>%
  dplyr::mutate(first2 = min(which(chloro == 0 | row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  dplyr::select(-first2)
length(df$ID)

# remove extreme chlorophyll point
# df <- head(df,-5)


sapply(split(df$lat,df$ID),length)
# drop the birds that have fewer than x relocations 
df <- df[!(as.numeric(df$ID) %in% which(table(df$ID)<2)),]
df <- droplevels(df)
length(df$ID)

# split up data by location - Ireland or Scotland 
df<-data.frame(df)
df <- df[df$location=="Scotland" , ] 
df<-droplevels(df)

# convert times from factors to POSIX
head(df)
df$date<-as.POSIXct(df$date, format= "%d/%m/%Y %H:%M", tz = "GMT")
head(df)

# prepare data with moveHMM 
# 8 day chloro
# stormData <- df[,c(4,5,8,13,17)]
# monthly chloro 
stormData <- df[,c(4,5,8,13,11)]

head(stormData)
stormData <- prepData(stormData,type="LL",coordNames=c("lon","lat"))
# plot(stormData,compact=T)

# add cosinor covariate based on hour of day
stormData$hour <- as.integer(strftime(stormData$date, format = "%H", tz="GMT"))
head(stormData)

# label 2 states
stateNames <- c("transiting", "foraging")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m1 <- list(step=c(10,5,1,2),angle=c(10,5)) # it's mean1,mean2,sd1,sd2 for step lengths

# fit model
m1 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m1,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)


# formula for transition probabilities with chloroMonth
formula <- ~ chloroMonth
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)


# formula for transition probabilities with chloroMonth quadratic 
formula <- ~chloroMonth+I(chloroMonth^2)
# initial parameters (obtained from nested model m1)
Par0_m3 <- getPar0(model=m2, formula=formula)
# fit model
m3 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0=Par0_m3$beta, stateNames = stateNames, formula=formula,retryFits = 1)



# label 3 states
stateNames <- c("transiting", "foraging", "resting")
# distributions for observation processes
dist = list(step = "gamma", angle = "vm")
# initial parameters
Par0_m4 <- list(step=c(10,5,1,1,2,1),angle=c(10,5,1))

# fit model
m4 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m4,
             estAngleMean = list(angle=FALSE), stateNames = stateNames, retryFits = 1)



# formula for transition probabilities with chloroMonth
formula <- ~ chloroMonth
# initial parameters (obtained from nested model m1)
Par0_m5 <- getPar0(model=m4, formula=formula)
# fit model
m5 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m5$Par,
             beta0=Par0_m5$beta, stateNames = stateNames, formula=formula,retryFits = 1)



# formula for transition probabilities with chloroMonth quadratic 
formula <- ~chloroMonth+I(chloroMonth^2) + ID
# initial parameters (obtained from nested model m1)
Par0_m6 <- getPar0(model=m5, formula=formula)
# fit model
m6 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m6$Par,
             beta0=Par0_m6$beta, stateNames = stateNames, formula=formula,retryFits = 1)




AIC(m1,m2,m3,m4,m5,m6)

# compute the pseudo-residuals
pr <- pseudoRes(m4)
# time series, qq-plots, and ACF of the pseudo-residuals
plotPR(m4,title("m4"))

# Saving R objects
saveRDS(m1, "m1Ire.rds")
saveRDS(m2, "m2Ire.rds")
saveRDS(m3, "m3Ire.rds")
saveRDS(m4, "m4Ire.rds")
saveRDS(m5, "m5Ire.rds")
saveRDS(m6, "m6Ire.rds")
saveRDS(m7, "m7Ire.rds")
saveRDS(m8, "m8Ire.rds")


plot(m4,covs = data.frame(hour = 4)) # 1-2


states <- viterbi(m4)
table(states)/length(states)
stormData$states <- states

boxplot(stormData$chloro~stormData$states)
boxplot(stormData$bath~stormData$states)

stormData %>% 
  as.data.frame() %>% 
  group_by(states, hour) %>% 
  dplyr::summarise(n = n()) %>% 
  ggplot(aes(hour, y = factor(states), fill = n)) + 
  geom_raster() +
  scale_fill_continuous('Count', low = 'white', high = 'darkgreen', guide = 'legend') +
  scale_x_continuous(expand = c(0, 0), breaks = c(1,seq(0, 25, 1), 24)) +
  ylab('State') + 
  xlab('Hour') +
  coord_fixed() +
  theme(legend.position = 'top',
        legend.key.height = unit(.2, 'cm'))