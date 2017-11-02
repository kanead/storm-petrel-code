###################################################################################
# Storm Petrel Movement Summary 
###################################################################################
# clean everything first
rm(list=ls())

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
mydata <- read.csv("movementSummary.csv",header = T , sep = ",")
head(mydata)
levels(mydata$stage)

# summary stats 
round (sapply(split(mydata$Max_range_.km. ,mydata$stage),mean),digits = 2)
round (sapply(split(mydata$Max_range_.km. ,mydata$stage),sd),digits = 2)

# run ANOVA
m1<-aov(Max_range_.km.~stage,data=mydata)
summary(m1)
boxplot(Max_range_.km.~stage,data=mydata) 
# examine residuals 
stdres = rstandard(m1)
qqnorm(stdres) 
qqline(stdres)

# apply kruskal-wallis if residuals non-normal 
m2<-kruskal.test(Max_range_.km.~stage,data=mydata)
m2

# Posthoc test 
library(FSA)
# stop R using scientific notation
options(scipen=999)
PT1 = dunnTest(Max_range_.km.~stage,data=mydata, method = "bh")
PT1

# ggplot 
p <- ggplot(mydata, aes(stage, Max_range_.km.))
p + geom_boxplot() + geom_jitter(width = 0.01)  + 
  theme_bw(base_size = 18) + scale_y_continuous(name = "max foraging range (km)") 
