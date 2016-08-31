library(deSolve)

SetparametersAndInitState <- function(krc, kcp, mp, D, k0, e1, e2, e3, f, fc, phi, initP){
    pd <- ifelse( D == 2 , pd2D, pd3D)
    a0 <- ifelse( D == 2 , a02D, a03D)
    h <- pv + 2*(D-1)*pd
    fm1 <- ifelse(fc == 1, "Gr", "Ac")
    fm2 <- ifelse(fc == 1, "Gr", ifelse(fc == 2, "Sw", "Ac"))
    fm3 <- ifelse(fc == 1, "Ac", ifelse(fc == 2, "Sw", "Ac"))
    krp <- krc * kcp
    f1 <- f(krc, D, pv, pd, fm1, phi)
    f2 <- f(krp, D, pv, pd, fm2, phi)
    f3 <- f(kcp, D, pv, pd, fm3, phi)

    ###intrapop parameters###
    r <- r0 * (krp * mp)**(b - 1)
    K <- k0 * (krp * mp)**(1 - b)
    q1 <- q0 * (kcp * mp)**(b - 1)
    q2 <- q0 * mp**(b - 1)
    ###interpop parameters###
    alfa1 <- a0 * (kcp * mp)**(h - 1) * f1
    alfa2 <- a0 * mp**(h - 1) * f2
    alfa3 <- a0 * mp**(h - 1) * f3
    ###Equilibrium abundances , receptor system
    Req <- q1 / (e1 * alfa1)
    Ceq <- (r / alfa1) * (1 - Req / K)
    
    params <- c(r = r, K = K, q1 = q1, q2 = q2, alfa1 = alfa1, alfa2 = alfa2, alfa3 = alfa3, e1 = e1, e2 = e2, e3 = e3)

    InitState <- c(R = Req, C = Ceq, P = initP)
    
    return(list(parameters = params, state = InitState))
}

LVIGP <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
        dR <- R * ( r * ( 1 - R / K) - (alfa1 * C  + alfa2 * P) )
        dC <- C * ( e1 * alfa1 * R - alfa3 * P - q1 )
        dP <- P * ( e2 * alfa2 * R + e3 * alfa3 * C - q2 )

        list(c(dR, dC, dP))
    })
}

SimulateInvPCR <- function(krc, kcp, mp, D, k0, e1, e2, e3, f, fc, phi, initP, finalT, steps, Dynamics){
    times <- seq(0, finalT, finalT / steps)
    Param_States <- SetparametersAndInitState(krc, kcp, mp, D, k0, e1, e2, e3, f, fc, phi, initP)
    out <- ode(y = Param_States$state, times = times, func = Dynamics, parms = Param_States$parameters)
    return(out)
}


Prueba <-  SimulateInvPCR(DD$K_RC[220], DD$K_CP[220] , DD$m_P[220], 3, 30, e1, e2, e3, f, DD$fc[220], 2.0, 0.001, 100000000, 10000, LVIGP)

plot(Prueba)
    



