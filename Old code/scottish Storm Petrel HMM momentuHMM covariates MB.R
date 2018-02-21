#--------------------------------------------------------------------------------
# Storm Petrel Movement Code Scottish Data Monthly Chlorophyll - Mark Bolton 
#--------------------------------------------------------------------------------
# clean everything first
rm(list=ls())
# careful not to load moveHmm alongside momentuHMM
library(momentuHMM)
library(rgdal)
library(dplyr)

# set the directory to wherever you've saved the models 
setwd("")

# Load the model objects
m1 <- readRDS("m1.rds")
m2 <- readRDS("m2.rds")
m3 <- readRDS("m3.rds")
m4 <- readRDS("m4.rds")
m5 <- readRDS("m5.rds")
m6 <- readRDS("m6.rds")
m7 <- readRDS("m7.rds")
m8 <- readRDS("m8.rds")

# you can plot the models with the following code and select the hour you want to show 
plot(m8,covs = data.frame(hour = 17)) # 1-2
plot(m8,covs = data.frame(hour = 20)) # 1-2 and 2-3

