#--------------------------------------------------------------------------------
# IDENTIFY POINTS LOCATED OVER LAND
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())

#load required libraries
library(maptools)

setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation- 8 Day Chlorophyll")
df<-read.csv("subsetTracks30MinInterpolation-monthlyChloro.csv",header=T,sep=",")
head(df)


## One example of a SpatialPolygons object mapping Earth's land areas
data(wrld_simpl)

## Create a SpatialPoints object
xy <- df[,c(4,5)]
pts <- SpatialPointsDataFrame(coords = xy, data = df,
                               proj4string=CRS(proj4string(wrld_simpl)))


## Find which points fall over land
ii <- !is.na(over(pts, wrld_simpl)$FIPS)
ii
## Check that it worked
plot(wrld_simpl)
points(pts, col=1+ii, pch=16)

df["land"] <- ii
write.table(df, file = "subsetTracks30MinInterpolationLand.csv", row.names=F, sep=",")
