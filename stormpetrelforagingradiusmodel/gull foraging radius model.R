rm(list=ls())

Er <- 6.8 # Watts - other activites 
V <- 11.1 # metres per second - flight speed 
Ej <- 1.2 # Watts - energy expenditure of juveniles
qf <- 0.05 # proprtion of time flying
L <- 200000 # load size in joules
Ef <- 10 # watts - energy expenditure during flight

birdRadius <- function(qf,Er,Ef,Ej,V,L){
E = (1-qf)*Er + (qf*Ef) + Ej
R = (V * L * qf) / (2 * E)
return(R)
}

birdRadius(qf,Er,Ef,Ej,V,L)


birdRadiusExpanded <- function(qf,Er,Ef,Ej,V,L){
  E = (1-qf)*Er + (qf*Ef) + Ej
  R = (V * L * qf) / (2 * ((1-qf)*Er + (qf*Ef) + Ej))
  return(R)
}

birdRadiusExpanded(qf,Er,Ef,Ej,V,L)


# Wilson's storm-petrel are mainly planktivorous birds 
# with a long gut passage time (>12 h)