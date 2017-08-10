# clean everything first
rm(list=ls())

library(ggplot2)
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data\\subsetTracks30MinInterpolation- 8 Day Chlorophyll")
stateData <- read.table("ScottishtracksPlusViterbi.csv", header=T,sep=",")
# stateData <- read.table("IrishtracksPlusViterbi.csv", header=T,sep=",")
head(stateData)

ggplot(stateData, aes(x=hour, y=viterbi, colour=factor(viterbi))) + 
  geom_count() +   scale_size_continuous(limits = c(1,30)) +
  scale_y_continuous(breaks = c(1,2,3)) +
  scale_colour_manual(values=c("1" = "goldenrod3", "2" = "lightskyblue", "3" = "seagreen3")) +
  theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


