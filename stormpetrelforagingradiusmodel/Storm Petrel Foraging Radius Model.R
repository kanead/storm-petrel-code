# Storm Petrel Foraging Radius Model 
rm(list=ls()) # remove everything currently held in the R memory
graphics.off() # close all open graphics windows 

# Wilson's Storm Petrel
e<-4330 # energy in joules per gram
m<-8 # stomach capacity in grams
t<-86400 # time at sea in seconds
r<-212000 # foraging radius in metres
v<-4.9 # cross country speed in metres per second
p<-.4 # power delivered to chick


power <- function(e,m,t,r,v) {
  
 (2*e*m)/(t + (2*r)/v)
}

power(e,m,t,r,v)

radius<-function(e,m,v,p,t){
((e*m*v)/p) - ((t*v)/2)
}

radius(e,m,v,p,t)

# European Storm Petrel
eSP<-4330
mSP<-5 # Diet and foraging behaviour of the british storm petrel in the bay of biscay during summer
tSP<-86400
rSP<-336000
vSP<-4.9

power(eSP,mSP,tSP,rSP,vSP)

pSP<-0.4

radius<-function(e,m,v,p,t){
  ((e*m*v)/p) - ((t*v)/2)
}

radius(eSP,mSP,vSP,pSP,tSP)

# use the function for radius where only one value can vary so we can plot it
radiusPlot<-function(p){
  (((4330*5*4.9)/p) - ((86400*4.9)/2))/1000
}
curve(radiusPlot, from=0.1, to=0.6, xlab="power (W)", ylab="radius (km)")
abline(v = 0.269, lty = 2)
