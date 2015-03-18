from sympy import exp

def set_r(r0,mr,w,Er,Tr,k):             
    return r0*mr**(w-1)*exp(-Er/(k*Tr))
def set_K(k0,mr,w,Ek,Tr,k):             
    return k0*mr**(1-w)*exp(Ek/(k*Tr))
def set_q1(q10,mc,w,Eq1,Tc,k):
    return q10*mc**(w-1)*exp(-Eq1/(k*Tc))
def set_q2(q20,mp,w,Eq2,Tp,k):
    return q20*mp**(w-1)*exp(-Eq2/(k*Tp))

