

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
    return ((q0/chi1)**(1/exp))
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
    val <- ifelse(fm == "Ac",sqrt(1 + k**(2*pv))*c*Pi,ifelse(fm == "Sw", (k**pv)*c*Pi , ifelse(fm == "Gr", c*Pi , NA)))
    return (val)
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

## Construct parameter set, for an initial list of parameter values
ConstructParamSet <- function(params){
    for(i in seq(2, length(params))){
        params = CombParams(params,i)
    }
    return(params)
}
##
# Assuming the same number of rows for all the list elements before i,
# add a new level of combination for the values in column i. The idea is
# to expand vertically each row until i to encompass all possible values of i
##
CombParams <- function(params, i){
    n <- length(params[[i]])
    params[[i]] <- rep(params[[i]],length(params[[1]]))
    for(j in seq(1,i-1)){
        params[[j]] <- as.vector(sapply(params[[j]], rep, n))
    }
    return(params)        
}

## Parameters especifications
D <- c(2,3)
phi <- c(0.1, 0.5, 1.0)
kLims <- c(-5,8)
kRange <- kLims[2]-kLims[1]
krc <- 10**seq(-5,8,kRange/n)
fm <- c("Ac", "Gr", "Sw")
# set parameter set
params <- list(D = D, phi = phi, krc = krc, fm = fm)
params <- ConstructParamSet(params)
k0 <- ifelse(params$D == 2 , 0.1, 30)
params$k0 = k0
# get mc
mc <-  InvasibilityRC(params$krc,params$D,params$k0,e,f,params$fm,params$phi)
mr <- krc*mc
# set main data frame for plotting
M <- data.frame(MassC =  mc , MassR = mr, krc = params$krc , Dim = factor(params$D), Prod = factor(params$k0), phi = factor(params$phi), ForagingMode = factor(params$fm))
levels(M$Dim) <- c("2D","3D")
levels(M$ForagingMode) <-  c("Captura Activa", "Pastoreo", "Captura Pasiva")

# plot
library(ggplot2)

Plot <-  ggplot(M,aes(x= log10(krc), y = log10(MassC),color = phi))  + geom_line() + facet_grid(Dim~ForagingMode) + xlab(expression(log[10](k[RC]))) + ylab(expression(log[10](m[C]))) + labs(colour = expression(phi)) + annotate("text", x = 0 , y = 6 , label= "frac(dC,dt) > 0",color = "black", fontface = "italic" , parse = TRUE) + theme_bw()

Plot
ggsave("R-CInv.pdf",limitsize=FALSE, height = 8, width = 8)  
