#--------------------------------------------------------------------------------
# Storm Petrel Movement Code 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)
library(dplyr)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation-bath")
df<-read.csv("subsetTracks30MinInterpolation-bath.csv",header=T,sep=",")

# rename columns
names(df)[names(df) == 'ETOPO1.Elevation'] <- 'bath'
names(df)[names(df) == 'location.long'] <- 'lon'
names(df)[names(df) == 'location.lat'] <- 'lat'
names(df)[names(df) == 'tag.local.identifier'] <- 'ID'
length(df$ID)

# remove positive bathymetry values
df<-group_by(df, ID) %>%
  dplyr:: mutate(first2 = min(which(bath > 0 | row_number() == n()))) %>%
  filter(row_number() <= first2) %>%
  select(-first2)
length(df$ID)

sapply(split(df$lat,df$ID),length)
# drop the birds that have fewer than x relocations 
df <- df[!(as.numeric(df$ID) %in% which(table(df$ID)<2)),]
df <- droplevels(df)
length(df$ID)
df<-data.frame(df)

# convert times from factors to POSIX
head(df)
df$date<-as.POSIXct(df$date, format= "%d/%m/%Y %H:%M", tz = "GMT")
head(df)

# split up data by location - Ireland or Scotland 
df <- df[df$lon>-3 , ] 
df<-droplevels(df)
df<-data.frame(df)

# prepare data with moveHMM
stormData <- df[,c(4,5,8,11,12)]
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
formula <- ~ bath
# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)
# fit model
m2 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula,retryFits = 1)


# formula for transition probabilities with time 
formula <- ~ bath* cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m3 <- getPar0(model=m1, formula=formula)
# fit model
m3 <- fitHMM(data = stormData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
             beta0=Par0_m3$beta, stateNames = stateNames, formula=formula,retryFits = 1)
#plot(m3)
#plot(m3,covs = data.frame(bath=3))
#plot(m3,covs = data.frame(hour=8))

# plot(m3,covs = data.frame(hour=10))
#plot(m3,covs = data.frame(bath=2, hour = 9))

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ bath* cosinor(hour, period = 24),
                       sd = ~ bath* cosinor(hour, period = 24)),
           angle = list(concentration = ~ bath))
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
formula <- ~ bath
# initial parameters (obtained from nested model m1)
Par0_m6 <- getPar0(model=m5, formula=formula)
# fit model
m6 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m6$Par,
             beta0=Par0_m6$beta, stateNames = stateNames, formula=formula,retryFits = 1)


# formula for transition probabilities with time 
formula <- ~ bath* cosinor(hour, period = 24)
# initial parameters (obtained from nested model m1)
Par0_m7 <- getPar0(model=m5, formula=formula)
# fit model
m7 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m7$Par,
             beta0=Par0_m7$beta, stateNames = stateNames, formula=formula,retryFits = 1)

# plot(m6, plotCI = FALSE, covs = data.frame(hour=7))

# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ bath* cosinor(hour, period = 24),
                       sd = ~ bath* cosinor(hour, period = 24)),
           angle = list(concentration = ~ bath))
# initial parameters (obtained from nested model m6)
Par0_m8 <- getPar0(model=m6, formula=formula, DM=DM)

# fit model
m8 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m8$Par,
             beta0 = Par0_m8$beta, DM = DM, stateNames = stateNames,
             formula = formula)


# formula for transition probabilities without time 
formula <- ~ ns(bath, knots = seq(0.1, 0.9, by = 0.1))
# initial parameters (obtained from nested model m1)
Par0_m9 <- getPar0(model=m5, formula=formula)
# fit model
m9 <- fitHMM(data = stormData, nbStates = 3, dist = dist, Par0 = Par0_m9$Par,
             beta0=Par0_m9$beta, stateNames = stateNames, formula=formula,retryFits = 1)



# plot(m8, plotCI = TRUE, covs = data.frame(hour=10))
AIC(m1,m2,m3,m4,m5,m6,m7,m8)

# compute the pseudo-residuals
pr <- pseudoRes(m6)
# time series, qq-plots, and ACF of the pseudo-residuals
plotPR(m6,title("m6"))

plot(m6,covs = data.frame(hour=20), plotCI=F)

states <- viterbi(m6)
stormData$states <- states
boxplot(stormData$bath~stormData$states)
boxplot(stormData$hour~stormData$states)

library(ggplot2)
# Basic violin plot
stormData$states <- as.factor(stormData$states)
p <- ggplot(stormData, aes(x=states, y=hour, fill=states)) + 
  geom_violin(trim=TRUE, adjust = 1) + geom_boxplot(width=0.1, fill="white") +
  scale_x_discrete(labels = c('exploratory','foraging','resting'))
p + theme_set(theme_gray(base_size = 20))+
coord_flip() + scale_fill_manual(values=c("#DAA520", "#87CEFA", "#006400"))+
theme(legend.position="none") # Remove legend

#goldenrod3
#lightskyblue
#darkgreen

# comparison of states
kruskal.test(stormData$hour~stormData$states)
# post hoc test
library(FSA)
PT1 = dunnTest(stormData$hour~stormData$states, method = "bh")
PT1

# Saving R objects
saveRDS(m1, "m1.rds")
saveRDS(m2, "m2.rds")
saveRDS(m3, "m3.rds")
saveRDS(m4, "m4.rds")
saveRDS(m5, "m5.rds")
saveRDS(m6, "m6.rds")
saveRDS(m7, "m7.rds")
saveRDS(m8, "m8.rds")
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation-bath\\Scottish data")
# Load the object
m1 <- readRDS("m1.rds")
m2 <- readRDS("m2.rds")
m3 <- readRDS("m3.rds")
m4 <- readRDS("m4.rds")
m5 <- readRDS("m5.rds")
m6 <- readRDS("m6.rds")
m7 <- readRDS("m7.rds")
m8 <- readRDS("m8.rds")
