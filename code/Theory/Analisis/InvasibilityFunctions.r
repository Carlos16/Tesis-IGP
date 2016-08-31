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
# Function to calculate the necessary conditions for the invasion of the top
# predator P, over the subsystem formed by R-C, asumming it was at equilibrium a
# at the moment of the invasion. if NecInvPCR <=0, then the invasion is not possible.
##

NecInvPCR <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    pd <- ifelse(D==2,pd2D,pd3D)
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    a0 <- e2*q0 / e1
    a1 <- e3*r0
    h <- pv + 2 * pd *(D-1)
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    f3 <- f(kcp, D, pv, pd, fm3, phi)
    chi3 <- a0 * f2 / f1
    chi4 <- a1 * f3 / (f1 * krc**(1-b))
    result <-  ( chi3 + chi4 ) * kcp**(b-h) - q0
    return (result)
}


InvPCR <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    Dim <- ifelse(c(D==2,D==2), c(pd2D, a02D), c(pd3D, a03D))
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    pd <- Dim[1]
    a0 <- Dim[2]
    h <- pv + 2*(D-1)*pd
    w <- h + 1 - 2*b
    ###############
    c0 <- e2 * q0 / e1
    c1 <- e3 * r0
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    f3 <- f(kcp, D, pv, pd, fm3, phi)
    chi3 <- c0 * f2 * kcp**(b - h) / f1
    chi4 <- c1 * f3 * kcp**(b - h) / ( f1 * krc**(1-b) )
    gamma1 <- kcp**(-w) * q0 / ( e1 * k0 * a0 * f1 * krc**(1-b))
    ################
    cond <-  chi3 + chi4 - q0
    cond2 <-  chi3 - q0
    coeff <- ifelse( cond2 > 0  , chi4/(chi4 + cond2) , 1)
    boundary <- ifelse(cond > 0 , (coeff * gamma1)**(1/w), NA)
    return (boundary)
}


NecInvCPR <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    pd <- ifelse(D == 2, pd2D , pd3D)
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    c0 <- e2 * q0 / e1
    h <- pv + 2 * pd *(D-1)
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    chi3 <- (c0 * f2 * kcp**(b - h)) / f1 
    result <- q0 - chi3
    return (result)
}

InvCPR <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    pd <- ifelse(D == 2, pd2D, pd3D)
    a0 <- ifelse(D == 2, a02D, a03D)
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    h <- pv + 2*(D-1)*pd
    w <- h + 1 - 2*b
    c <-  e2/(e1 * e3)
    ###############
    c0 <- e2 * q0 / e1
    c1 <- e3 * r0
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    f3 <- f(kcp, D, pv, pd, fm3, phi)
    chi3 <- (c0 * f2 * kcp**(b - h)) / f1
    chi4 <- (c1 * f3 * kcp**(b - h)) / ( f1 * krc**(1-b) )
    gamma2 <- q0 / ( e2 * k0 * a0 * f2 * krp**(1-b))
    ################
    cond1 <- chi3 + c * chi4 - q0
    cond2 <- chi3 - q0

    boundary2 <- ifelse(cond2 < 0 , (c * chi4 * gamma2 /cond1)**(1/w) , NA)
    boundary <-  ifelse(cond1 <=0 , NA, boundary2)
    return (boundary)
}

SufEstabCoex <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    pd <- ifelse(D == 2, pd2D, pd3D)
    a0 <- ifelse(D == 2, a02D, a03D)
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    h <- pv + 2*(D-1)*pd
    w <- h + 1 - 2*b
    c <-  e2/(e1 * e3)
    ###############
    c1 <- e3 * r0
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    f3 <- f(kcp, D, pv, pd, fm3, phi)
    chi4 <- (c1 * f3 * kcp**(b - h)) / ( f1 * krc**(1-b) )
    chi2 <-  e2 * k0 * a0 * f2 * krp**(1-b)
    result <- ((c/(c-1)) * chi4 / chi2)**(1/w)
    return(result)
        
}

CondCPR <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    Dim <- ifelse(c(D==2,D==2), c(pd2D, a02D), c(pd3D, a03D))
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    pd <- Dim[1]
    a0 <- Dim[2]
    h <- pv + 2*(D-1)*pd
    w <- h + 1 - 2*b
    c <-  e2/(e1 * e3)
    ###############
    c0 <- e2 * q0 / e1
    c1 <- e3 * r0
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    f3 <- f(kcp, D, pv, pd, fm3, phi)
    chi3 <- c0 * f2 * kcp**(b - h) / f1
    chi4 <- c1 * f3 * kcp**(b - h) / ( f1 * krc**(1-b) )
    ################
    cond <-  chi3 + c * chi4 - q0
    result <-  ifelse(cond <= 0 , TRUE, FALSE)
}


InvPR <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
    Dim <- ifelse(c(D==2,D==2), c(pd2D, a02D), c(pd3D, a03D))
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    pd <- Dim[1]
    a0 <- Dim[2]
    h <- pv + 2*(D-1)*pd
    w <- h + 1 - 2*b
    ###############
    krp <- krc * kcp
    f2 <- f(krp, D, pv, pd, fm2, phi)
    gamma2 <- q0 / ( e2 * k0 * a0 * f2 * krp**(1-b))
    ################
    boundary <- gamma2**(1/w)
    return (boundary)
}



getZone <- function(mass, lowbound, highbound, cond){
    r1 <- mass > lowbound
    return ( ifelse(cond, r1 , r1 & mass < highbound))
}


CondInvCoex <- function(krc, kcp, mp, D, k0, e1, e2, e3, f, fc, phi1, phi2){
    pd <- ifelse(D == 2, pd2D, pd3D)
    a0 <- ifelse(D == 2, a02D, a03D)
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    h <- pv + 2*(D-1)*pd
    w <- h + 1 - 2*b
    c <-  e2/(e1 * e3)
    ###############
    c0 <- e2 * q0 / e1
    c1 <- e3 * r0
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi1)
    f2 <- f(krp, D, pv, pd, fm2, phi2)
    f3 <- f(kcp, D, pv, pd, fm3, phi2)
    chi1 <- e1 * k0 * a0 * f1 * krc**(1-b)
    chi3 <- (c0 * f2 * kcp**(b - h)) / f1
    chi4 <- (c1 * f3 * kcp**(b - h)) / ( f1 * krc**(1-b) )
    chi2 <-  e2 * k0 * a0 * f2 * krp**(1-b)
    gamma2 <- q0 / ( e2 * k0 * a0 * f2 * krp**(1-b))
    gamma1 <- kcp**(-w) *q0 / chi1
    condPCR <- chi3 + chi4 - q0
    condCPR <- chi3 + c*chi4 - q0
    NCondCPR <- q0 - chi3
    gamma3 <- chi4/condPCR * gamma1
    gamma4 <- (c*chi4/condCPR) * gamma2
    #Am <- ((c/(c-1)) * chi4 / chi2)**(1/w)
    
    ################
    mu1 <- mp > gamma1**(1/w)
    mu2 <- mp > gamma2**(1/w)
    om1 <- condPCR > 0
    mgzeta3 <-  mp > gamma3**(1/w)
    mu3 <- mgzeta3 & om1
    om2 <- condCPR < 0
    mszeta4 <- mp < gamma4**(1/w)
    mu4 <- om2 | mszeta4
    #Est <- mp < Am
    InvPCR <-  mu1 & mu3
    InvCPR <-  mu2 & mu4
    Coex <- mu3 & mu4
    #UnEstCoex <- !Est & !mu3 & !mu4
    MutInv <- InvPCR & InvCPR
    FulCond <- list(mu1, mu2, mu3, mu4, InvPCR, InvCPR, Coex, MutInv)
    return(FulCond)
}


GetProp <- function(Data, D, k0, e1, e2, e3, f, phi1, phi2){    
    T <- with(Data,CondInvCoex(K_RC, K_CP, m_P, D, k0, e1, e2, e3, f, fc, phi1, phi2))
    n <- nrow(Data)
    Num <- as.numeric(lapply(T, sum , na.rm = T))
    return(Num/n)
}

getPropk <- function(k){
    return (GetProp(Data, D, k, e1, e2, e3, f, phi1, phi2))
}
