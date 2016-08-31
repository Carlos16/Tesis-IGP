MTPCond <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
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
    ################
    cond <- chi3/c + chi4 - q0
    cond2 <- q0 - chi3
        
    return (ifelse(cond2 > 0 , cond , NA))
}
