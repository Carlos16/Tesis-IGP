CoexDiff <- function(krc, kcp, D, k0, e1, e2, e3, f, fc, phi){
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
    chi1 <- e1 * k0 * a0 * f1 * krc**(1-b)
    chi3 <- (c0 * f2 * kcp**(b - h)) / f1
    chi4 <- (c1 * f3 * kcp**(b - h)) / ( f1 * krc**(1-b) )
    chi2 <-  e2 * k0 * a0 * f2 * krp**(1-b)
    gamma2 <- q0 / ( e2 * k0 * a0 * f2 * krp**(1-b))
    gamma1 <- kcp**(-w) *q0 / chi1
    condPCR <- chi3 + chi4 - q0
    condCPR <- chi3 + c*chi4 - q0
    gamma3 <- ifelse(condPCR > 0 , log10(chi4 / condPCR * gamma1), NA)
    gamma4 <- ifelse(condCPR > 0 , log10(c * chi4 / condCPR * gamma2), NA )
    diff <- gamma4 - gamma3
    return (ifelse( diff > 0 , diff/w , NA))
}


