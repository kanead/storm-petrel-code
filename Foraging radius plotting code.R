Er <- 0.68#  Watts - energetic cost of other activites 
V <- 4.9 #  metres per second - flight speed 
Ej <- 0.95 # Watts - energy expenditure of juveniles
qf <- 0.95# 0.05 # proprtion of time flying
L <-  21650 #  load size in joules * energy density of food
Ef <- 1.04 # watts - energy expenditure during flight
Me <- 6 # number of meals the bird eats 

petrelRadius <- function(qf,Er,Ef,Ej,V,L,Me){
  R = (V * (L*Me) * qf) / (2 * ((1-qf)*Er + (qf*Ef) + Ej))
  return(R/1000) # return in km 
}

petrelRadius(qf,Er,Ef,Ej,V,L,Me)


# start
petrelRadius <- function(Me){
  R = (4.9 * 21650 * Me * 0.95) / (2 * ((1-0.95)*0.68 + (0.95*1.04) + 0.95))
  return(R/1000) # return in km 
}
curve(petrelRadius, from=1, to=6, xlab="Feeding Events", ylab="Radius (km)", ylim=c(0,300),lty=2)



petrelRadius2 <- function(Me){
  R = (4.9 * 21650 * Me * 0.95) / (2 * ((1-0.95)*0.68 + (0.95*1.04) + 0))
  return(R/1000) # return in km 
}
# curve(petrelRadius, from=1, to=6, xlab="Feeding Events", ylab="radius (km)")

plot (petrelRadius2, 1, 6, add=TRUE)
abline(h=210, col="gold2", lty=1)
abline(h=161, col="red4", lty=1)
abline(h=228, col="chartreuse4", lty=1)

# 209.78 incubating 
# 161.15 brooding
# 228.00 post brooding 


