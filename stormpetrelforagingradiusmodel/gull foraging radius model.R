# Foraging Radius Model 


rm(list=ls())
# Bird 908
# 66.9 hours flying 
# Avg speed = 4.6m/s
# max speed = 11.2 m/s
# 50% of time transiting

Er <- 0.6# 6.8 # Watts - other activites 
V <- 4.9 # 11.1 # metres per second - flight speed 
Ej <- 0.6# 1.2 # Watts - energy expenditure of juveniles
qf <- 0.75# 0.05 # proprtion of time flying
L <-  21650 # 200000 # load size in joules
Ef <- 1.42 # 10 # watts - energy expenditure during flight

birdRadius <- function(qf,Er,Ef,Ej,V,L){
E = (1-qf)*Er + (qf*Ef) + Ej
R = (V * L * qf) / (2 * E)
return(R/1000) # return in km 
}

birdRadius(qf,Er,Ef,Ej,V,L)


birdRadiusExpanded <- function(qf,Er,Ef,Ej,V,L){
  R = (V * L * qf) / (2 * ((1-qf)*Er + (qf*Ef) + Ej))
  return(R/1000) # return in km 
}

birdRadiusExpanded(qf,Er,Ef,Ej,V,L)


# Wilson's storm-petrel are mainly planktivorous birds 
# with a long gut passage time (>12 h)

# References 
# -  How to Live in Colonies: Foraging Range and Patterns of Density around a Colony of 
# Black-Headed Gulls Larus ridibundus in Relation to the Gulls' Energy Budget
# -  Daily energy expenditure by adult leach's storm petrels during the nesting cycle 