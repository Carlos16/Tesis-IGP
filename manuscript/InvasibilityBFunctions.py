import numpy as np
import matplotlib.pyplot as plt
import matplotlib as m


def g(k_,pv,pd,T1,T2,E1,E2,D,v01,v02,fm,thermy1,thermy2,k):
    if thermy1=="Ecto" and thermy2=="Ecto":
        if fm=="Active":
            delta=np.exp(-(1./k)*(E1/T1 - E2/T2))
            return v02*np.exp(-E2/(k*T2))*(1+(v01/v02)*k_**(2*pv)*delta**2)**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)*np.exp(-E1/(k*T1))
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)*np.exp(-E2/(k*T2))
    elif thermy1=="Endo" and thermy2=="Ecto":
        if fm=="Active":
            return v02*(1+(v01/v02)*k_**(2*pv)*np.exp(-2*E1/(k*T1)))**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)*np.exp(-E1/(k*T1))
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)
    elif thermy1=="Ecto" and thermy2=="Endo":
        if fm=="Active":
            return v02*(1+(v01/v02)*k_**(2*pv)*np.exp(2*E2/(k*T2)))**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)*np.exp(-E2/(k*T2))

def f(k_,form,a,b,c):
    if form == 1:
        return a
    elif form == 2:
        return a/(1+k_**b)
    elif form == 3 :
        return a*exp(-(log(k_)-b)**2/(2*c))




def InvP3(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KRC*KCP
    f2 = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    n0 = e2*k0*alfa0*f2*KRP**(1-w)
    h = pv + 2*(D-1)*pd 
    v = h + 1 - 2*w
    return (1./v)*np.log10(q20/n0)

def InvC2(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    h = pv + 2*(D-1)*pd 
    KRP = KRC*KCP
    f1 = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    x0 = e1*k0*alfa0*f1*KCP**(h-1)*KRP**(1-w)
    x1 = q10*KCP**(w-1)
    v = h + 1 - 2*w
    return (1./v)*np.log10(x1/x0)

def InvP5(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    g_C = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    g_R = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    g_P = g(KCP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpc,'Ecto','Ecto',1.)*f(KCP,2,1,b,0)
    h = pv+2*(D-1)*pd
    R0 = q20/(e2*alfa0*g_R)
    C0 = e1*alfa0*g_C*KCP**(h-1)*R0
    P_1 = r0*KRP**(2*(w-1))/(e2*k0*(alfa0*g_R)**2)
    P_2 = e2*alfa0*g_R*k0*KRP**(1-w)
    C1 = alfa0*g_P*P_1*P_2 
    C2 = alfa0*g_P*P_1*q20
    C3 = q20*KCP**(w-1)
    
    #print ((C1 + C3 - C0)/C2)
    return (1./(2*w -(h +1)))*np.log10(( C3 + C1 - C0)/C2)


def Dbound(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    f1 = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    f2 = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    f3 = g(KCP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpc,'Ecto','Ecto',1.)*f(KCP,2,1,b,0)
    h = pv+2*(D-1)*pd
    D_e = e2 - e1*e3
    u1 = k0*KRP**(2*(1-w))*alfa0**2*f1*f2*KCP**(h-1)*D_e
    u0 =  r0*e3*alfa0*f3
    
    #print ((C1 + C3 - C0)/C2)
    return (1./(h +1 - 2*w))*np.log10(u0/u1)

def InvP5B(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    g_C = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    g_R = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    g_P = g(KCP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpc,'Ecto','Ecto',1.)*f(KCP,2,1,b,0)
    h = pv+2*(D-1)*pd
    R0 = q20/(e2*alfa0*g_R)
    C0 = e1*alfa0*g_C*KCP**(h-1)*R0
    P_1 = r0*KRP**(2*(w-1))/(e2*k0*(alfa0*g_R)**2)
    P_2 = e2*alfa0*g_R*k0*KRP**(1-w)
    C1 = alfa0*g_P*P_1*P_2 
    C2 = alfa0*g_P*P_1*q20
    C3 = q20*KCP**(w-1)
    
    #print ((C1 + C3 - C0)/C2)
    return  C3 + C1 - C0


def InvP4(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    f1 = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    f2 = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    f3 = g(KCP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpc,'Ecto','Ecto',1.)*f(KCP,2,1,b,0)
    h = pv+2*(D-1)*pd
    t0 = q10*KCP**(w-h)/(e1*alfa0*f1)
    t1 = r0*KRP**(w-1)/(alfa0*f1*KCP**(h-1))
    t2 = t0/(k0*KRP**(1-w))
    t3 = e2*alfa0*f2*t0
    t4 = e3*alfa0*f3*t1
    #print ((C1 + C3 - C0)/C2)
    return (1./(h +1 - 2*w))*np.log10((t4*t2)/(t3 + t4 - q20))




def InvP4B(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    f1 = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    f2 = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    f3 = g(KCP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpc,'Ecto','Ecto',1.)*f(KCP,2,1,b,0)
    h = pv+2*(D-1)*pd
    t0 = q10*KCP**(w-h)/(e1*alfa0*f1)
    t1 = r0*KRP**(w-1)/(alfa0*f1*KCP**(h-1))
    t2 = t0/(k0*KRP**(1-w))
    t3 = e2*alfa0*f2*t0
    t4 = e3*alfa0*f3*t1
    #print ((C1 + C3 - C0)/C2)
    return t3 + t4 - q20


def InvP4B2(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    f1 = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    f2 = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    f3 = g(KCP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpc,'Ecto','Ecto',1.)*f(KCP,2,1,b,0)
    h = pv+2*(D-1)*pd
    t0 = q10*KCP**(w-h)/(e1*alfa0*f1)
    t1 = r0*KRP**(w-1)/(alfa0*f1*KCP**(h-1))
    t2 = t0/(k0*KRP**(1-w))
    t3 = e2*alfa0*f2*t0
    #print ((C1 + C3 - C0)/C2)
    return t3 - q20


def Rrule(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    KRP = KCP*KRC
    f1 = g(KRC,pv,pd,1.,1.,0.,0.,D,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b,0)
    f2 = g(KRP,pv,pd,1.,1.,0.,0.,D,1.,1.,fpr,'Ecto','Ecto',1.)*f(KRP,2,1,b,0)
    h = pv+2*(D-1)*pd
    C1 = q20/(e2*alfa0)
    C2 = q10/(e1*alfa0)
    
    
    #print ((C1 + C3 - C0)/C2)
    return f1 - (C2/C1 )* f2*KCP**(w-h)
