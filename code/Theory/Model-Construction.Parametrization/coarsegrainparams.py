from sympy import exp
b_k= 8.6173324*10**(-5) #boltzmann constant in ev


"""metabolic and biomechanic parameterization of the
 search rate taking into consideration foraging strategies"""

def set_th(t_h0,m,w,E,k,T):
    r"""
    
    """
    return t_h0*m**(-w)*exp(-E/(k*T))

def set_alfa0(d0,D):
    if D==2:
        return 2*np.pi*d0
    elif D==3:
        return 4*np.pi*d0**2

def alfa(m2,alfa0,pv,pd,D,g,f):
    return alfa0*m2**(pv+2*(D-1)*pd)*g*f

def g(k_,pv,pd,T1,T2,E1,E2,D,v01,v02,fm,thermy1,thermy2,k):
    if thermy1=="Ecto" and thermy2=="Ecto":
        if fm=="Active":
            delta=exp(-(1./k)*(E1/T1 - E2/T2))
            return v02*exp(-E2/(k*T2))*(1+(v01/v02)*k_**(2*pv)*delta**2)**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)*exp(-E1/(k*T1))
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)*exp(-E2/(k*T2))
    elif thermy1=="Endo" and thermy2=="Ecto":
        if fm=="Active":
            return v02*(1+(v01/v02)*k_**(2*pv)*exp(-2*E1/(k*T1)))**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)*exp(-E1/(k*T1))
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)
    elif thermy1=="Ecto" and thermy2=="Endo":
        if fm=="Active":
            return v02*(1+(v01/v02)*k_**(2*pv)*exp(2*E2/(k*T2)))**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)*exp(-E2/(k*T2))

def f(k_,form,a,b,c):
    if form == 1:
        return a
    elif form == 2:
        return a/(1+k_**b)
    elif form == 3 :
        return a*exp(-(log(k_)-b)**2/(2*c))

def set_alfa(m2,alfa0,k_,pv,pd,T1,T2,E1,E2,D,v01,v02,g,alfa,fm,thermy1,thermy2,k,a,b,c,form):
    g_= g(k_,pv,pd,T1,T2,E1,E2,D,v01,v02,fm,thermy1,thermy2,k)
    f_ = f(k_,form,a,b,c)
    alfa_=alfa(m2,alfa0,pv,pd,D,g_,f_)
    return alfa_


##Intra population parameters
def set_r(r0,mr,w,Er,Tr,k):             
    return r0*mr**(w-1)*exp(-Er/(k*Tr))
def set_K(k0,mr,w,Ek,Tr,k):             
    return k0*mr**(1-w)*exp(Ek/(k*Tr))
def set_q1(q10,mc,w,Eq1,Tc,k):
    return q10*mc**(w-1)*exp(-Eq1/(k*Tc))
def set_q2(q20,mp,w,Eq2,Tp,k):
    return q20*mp**(w-1)*exp(-Eq2/(k*Tp))



