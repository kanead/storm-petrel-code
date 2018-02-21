# MarMap Package
# http://www.molecularecologist.com/2015/07/marmap/

# clean everything first
rm(list=ls())

# Load package
library(marmap)
library(ggplot2)

# load tracking data
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
#setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\Storm Petrels\\Tracking Data")
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

# select Scottish data only 
scottishdata <- data[data$location=="scotland" , ] 
scottishdata<-droplevels(scottishdata)

# rename it 
track1<-scottishdata
# match it to example data
myvars <- c("lon", "lat","ID")
track1 <- track1[myvars]
head(track1)

# select the Irish data only 
irishdata <- data[data$location=="ireland" , ] 
irishdata<-droplevels(irishdata)
irishdata<-irishdata[irishdata$lat < 54.5, ]

# rename it 
track2<-irishdata
# match it to example data
myvars <- c("lon", "lat","ID")
track2 <- track2[myvars]
head(track2)

library(mapproj)
# Download bathymetric data and save on disk - Irish 
batIre <- getNOAA.bathy(-16,-5.5, 51, 56, res = 1, keep = TRUE)

# Download bathymetric data and save on disk - Scottish 
batScot <- getNOAA.bathy(-8,6, 56, 62, res = 1, keep = TRUE)

# Get depth profile along both tracks and remove positive depths (since random fake values can be on land)
# path.profile() gets depth value for all cells of the bathymetric grid below the gps tracks
path1 <- path.profile(track1[,-3], batScot) ; path1 <- path1[-path1[,4]>0,]

# Get depth values for each gps tracking point
# get.depth() retrieve depth values only for gps tracking points
depth1 <- get.depth(batScot, track1$lon, track1$lat, locator = FALSE, distance = TRUE) 

# Add depth values to tracks 1 and 2 and remove positive depths
track1$depth <- depth1$depth ; track1 <- track1[-track1$depth > 0,]

# Plot
# layout(matrix(c(1, 2, 1, 3, 4, 4), ncol = 2, byrow = TRUE), height = c(1, 1, 1))

# Select colours
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

## Bathymetric map with gps tracks - Scottish 
plot(batScot, land = TRUE, image = TRUE, lwd = 0.2, bpal = list(c(min(batScot,na.rm=TRUE), 0, blues), c(0, max(batScot, na.rm=TRUE), greys)), ylim=c(56,62))
plot(batScot, deep = 0, shallow = 0, step = 0, lwd = 0.8, add = TRUE)

# Scottish Tracks 
color_pallete_function <- colorRampPalette(
  colors = c("red", "orange", "blue"),
  space = "Lab" # Option used when colors do not represent a quantitative scale
)

autoplot(batScot, geom=c("r", "c")) +
  scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen")
# p1$plot + geom_line(data = track1, mapping = aes_string(x='lon', y='lat'), inherit.aes = F)

last_plot() + geom_point(aes(x=lon, y=lat, colour=factor(ID)), data=track1, alpha=0.5)

last_plot() + theme(legend.position="none")

last_plot() + labs(x = "longitude", y = "latitude")

last_plot() + scale_fill_etopo()

# set aesthetics
#autoplot(batScot, geom=c("r", "c"), colour="white", size=0.1)

# topographical colour scale, see ?scale_fill_etopo
#autoplot(batScot, geom=c("r", "c"), colour="white", size=0.1) + scale_fill_etopo()

##############################################################################################################
## Bathymetric map with gps tracks - Irish
##############################################################################################################
plot(batIre, land = TRUE, image = TRUE, lwd = 0.2, bpal = list(c(min(batIre,na.rm=TRUE), 0, blues), c(0, max(batIre, na.rm=TRUE), greys)), ylim=c(51,56))
plot(batIre, deep = 0, shallow = 0, step = 0, lwd = 0.8, add = TRUE)


# Irish Tracks 
color_pallete_function <- colorRampPalette(
  colors = c("red", "orange", "blue"),
  space = "Lab" # Option used when colors do not represent a quantitative scale
)

autoplot(batIre, geom=c("r", "c")) +
  scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen")
# p1$plot + geom_line(data = track1, mapping = aes_string(x='lon', y='lat'), inherit.aes = F)

last_plot() + geom_point(aes(x=lon, y=lat, colour="factor(states)"), data=tracks2, alpha=0.5)
last_plot() + geom_point(aes(x=lon, y=lat, colour="#00BA80"), data=tracks2, alpha=0.5)


last_plot() + geom_point(aes(x=lon, y=lat), data=track2, alpha=0.5)

last_plot() + theme(legend.position="none")

last_plot() + labs(x = "longitude", y = "latitude")

last_plot() + scale_fill_etopo()

last_plot() + theme_black() 

head(tracks)
tracks2 <- tracks[tracks$states=="2",]
tracks2 <- droplevels(tracks2)
theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}


