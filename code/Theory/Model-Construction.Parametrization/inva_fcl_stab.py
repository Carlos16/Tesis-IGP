#Calculating food chain length within the stable coexistence region
#LV case

def setI_R_LV(R,a2,e2,m_P):
    ha2 = a2/m_P
    return e2*ha2*R

def setI_C_LV(C,a3,e3,m_P):
    ha3 = a3/m_P
    return e3*ha3*C

#LV case at equilibrium
def set_MTP_C(Req,a2,m_P,q2,e2):
    return 3. - (a2/m_P)*e2*Req/q2


#RM case

def setI_R_RM(R,e2,a2,m_P,t_hp):
    ha2 = a2/m_P
    return e2*ha2*R/(1+t_hp*ha2*R)

def setI_C_RM(C,e3,a3,m_P,t_hp):
    ha3 = a3/m_P
    return e3*ha3*C/(1+t_hp*ha3*C)
#RM case at equilibrium
def setMPTEqRM(I_R,q2):
    return 3  - I_R/q2

#MTP in general for both cases
def setMTP(I_R,I_C):
    return (2*I_R + 3*I_C)/(I_R+I_C)

#Invasibility criterions

##C to R , scenario 2


###Lotka-Volterra

def set_I_C_s2(e1,a1,K,q1):
    return e1*a1*K - q1
###Rosenzweig-MacArthur

def set_I_C_s2RM(e1,a1,K,q1,hC0,q10):
    return a1*K*(e1 - q10*hC0) - q1

##P to R , scenario 3

###Lotka-Volterra
def set_I_P_s3(e2,a2,K,q2):
    return e2*a2*K - q2

###Rosenzweig-MacArthur
def set_I_P_s3RM(e2,a2,K,q2,hP0,q20):
    return a2*K*(e2-q20*hP0) - q2

##P to C-R , scenario 4

###Lotka-Volterra
def set_I_P_s4(e2,e3,a2,a3,q2,R,C):
    return e2*a2*R + e3*a3*C - q2

###Rosenzweig-MacArthur
def set_I_P_s4RM(e2,e3,a2,a3,q2,R,C,hP0,q20):   
    return a2*R*(e2 - q20*hP0) + a3*C*(e3 - q20*hP0)  - q2

##C to P-R , scenario 5

###Lotka-Volterra
def set_I_C_s5(e1,a1,a3,R,P,q1):
    return e1*a1*R - a3*P - q1

###Rosenzweig-MacArthur
def set_I_C_s5RM(e1,e2,a1,a3,m_C,R,P,q1,t_hc,q10,q20,hP0,hC0):
    return (e2/(e2 - q20*hP0)) *(a1*R*(e1 - q10*hC0) - q1) - a3*P*(1+t_hc*a1*m_C*R)



# Functions for calculating Req and Ceq in scenario 2 , Req and Peq in scenario 3 are similar

##Lotka-Volterra
def set_R_C_eq_sLV(r,K,q1,a1,e1):
    Req = q1/(e1*a1)
    Ceq = (r/a1)*(1- q1/(e1*a1*K))
    return Req, Ceq

##Rosenzweig-MacArthur
def set_R_C_eq_sRM(r,K,q1,q10,a1,e1,hC0):
    Req = q1 /(a1*(e1 - q10*hC0))
    Ceq = (e1*Req*r/q1)*(1-Req/K)
    return Req,Ceq

#Qualitative stability analysis using Hurwtiz determinants

##auxiliary functions
def set_D(K,a1,a2,a3,e1,e2,e3,m_C,r):
    return K*a1*a2*(e1*e3-e2) + a3*e3*m_C*r

def set_d1(r,Req,K):
    return  r*Req/K

def set_d2(e1,e2,e3,a1,a2,a3,m_C,m_P,Ceq,Req,Peq):
    return (e1*(a1*m_P)**2*Ceq*Req + (Peq*m_C**2)*(e2*a2**2*Req + e3*a3**2*Ceq))/(m_C*m_P)**2

def set_d3(D,a3,Ceq,Req,Peq,K,m_C,m_P):
    return a3*Req*Ceq*Peq*D/(K*m_C*m_P**2)

##Hurwtiz criterion
def set_hdet2(d1,d2,d3):
    return d1*d2-d3


