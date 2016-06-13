

## Biomecanic constants

pv = 0.26
pd2D = 0.21
pd3D = 0.2
a02D = 10**(-3.08)
a03D = 10**(-1.77)

## Metabolic constant

b = 0.75
q0 = 4.15e-8

##
# Function to calculate the necessary mass for the consumer to invade a system in which is resource
# is at carrying capacity K. It uses the relationships derived in Pawar et.al 2012
# D : Dimension of the space in which the interaction takes place
# phi: measures the decay in the probability of capture following an encounter with respect to increases in the size ratio
# krc: size ratio resource / consumer
# f: function that captures the influence of the size ratio over the search rate of the consumer
# fm: foraging mode of the predator
# e: conversion efficiency  biomass resource \to biomass predator
# k0: basal productivity of the habitat, take note that changes in k0 \to ck0, moves the graph of (krc,Inv(krc)) by :
# -log10(c)*1/(hr + 1 - 2 *b)
##
InvasibilityRC <- function(krc, D, k0, e, f, fm, phi){
    Dim <- ifelse(c(D==2,D==2),c(pd2D, a02D), c(pd3D, a03D))
    pd <- Dim[1]
    a0 <- Dim[2]
    hr = pv + 2*pd*(D-1)
    chi1 = e*k0*a0*f(krc, D, pv, pd, fm, phi)*krc**(1-b)
    exp = hr + 1 - 2*b
    return (q0/chi1)**(1/exp)
}


##
#Function that computes the influence of the size ratios over the search rate of the consumer
#fm: foraging mode, available modes : "Active", "Passive" and "Sit-Wait"
#D: Dimension of the habitat.
#pd,pv: biomechanical constants, pd is influences by D
#phi: as defined in InvasibilityRC
##
f <- function(k, D, pv, pd, fm, phi){
    Pi = 1 / (1 + k**phi)
    c = k**((D-1)*pd)
    if(fm == "Ac"){
        return (sqrt(1 + k**(2*pv))*c*Pi)
    }else if(fm == "Gr"){
        return ((k**pv)*c*Pi)
    }else if (fm == "Sw"){
        return (c*Pi)
    }else{
        print("NOT VALID FORAGING MODE")
    }
}

##
##Auxiliary function to set the initial D and k0 vectors                                     
##
assingD <- function(i,n){
    if( i <= n)
        return (2)
    else
        return (3)
}

e = 0.3
n = 700
k0 = 300
D = 3
fm = "Ac"
phi = 5

## Construct parameter set, modify it for foraging mode instead of k0

k0Base <- c(0.01, 0.1, 1, 3, 30, 300)
n2 <- length(k0Base)
D <- sapply(seq(1:n2),assingD,n2/2)

phiBase <- c(1,3,5)
n3 <- length(phiBase)
phi <- rep(phiBase,n2)
D <- as.vector(sapply(D,rep,n3))
k0 <- as.vector(sapply(k0Base,rep,n3))

kLims <- c(-5,8)
kRange <- kLims[2]-kLims[1]

krcBase <- 10**seq(-5,8,kRange/n)
krc <-  rep(krcBase,n2*n3)

phi <- as.vector(sapply(phi,rep,n+1))
D <- as.vector(sapply(D,rep,n+1))
k0 <- as.vector(sapply(k0,rep,n+1))

mc <-  InvasibilityRC(krc,D,k0,e,f,fm,phi)
mr <- krc*mc

M <- data.frame(MassC =  mc , MassR = mr, krc = krc , Dim = factor(D), Prod = factor(k0), phi = factor(phi))

library(ggplot2)

Plot <-  ggplot(M,aes(x= log10(krc), y = log10(MassC),color = Prod))  + geom_line() + facet_grid(Dim~phi)
Plot

