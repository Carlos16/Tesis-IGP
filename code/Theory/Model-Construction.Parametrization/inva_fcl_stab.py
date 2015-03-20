#Calculating food chain length within the stable coexistence region
#LV case

def setI_R_LV(R,a2,e2,m_P):
    ha2 = a2/m_P
    return e2*ha2*R

def setI_C_LV(C,a3,e3,m_P):
    ha3 = a3/m_P
    return e3*ha3*C

#LV case at equilibrium
def set_MTP_C(Req,a2,m_P,q_2,e2):
    return 3. - (a2/m_P)*e2*Req/q_2


#RM case

def setI_R_RM(R,e2,a2,m_P,t_hp):
    ha2 = a2/m_P
    return e2*ha2*R/(1+t_hp*ha2*R)

def setI_C_RM(C,e3,a3,m_P,t_hp):
    ha3 = a3/m_P
    return e3*ha3*C/(1+t_hp*ha3*C)
#RM case at equilibrium
def setMPTEqRM(I_R,q_2):
    return 3  - I_R/q_2

#MTP in general for both cases
def setMTP(I_R,I_C):
    return (2*I_R + 3*I_C)/(I_R+I_C)

#Invasibility criterions LV model

##C to R , scenario 2
def set_I_C_s2(e1,alfa_1,m_C,K,q_1):
    return e1*(alfa_1/m_C)*K - q_1

##P to R , scenario 3
def set_I_P_s3(e_2,alfa_2,m_P,K,q_2):
    return e_2*(alfa_2/m_P)*K - q_2

##P to C-R , scenario 4
def set_I_P_s4(e_2,e_3,alfa_2,alfa_3,m_P,q_2,R,C):
    return e_2*alfa_2*R + e_3*alfa_3*C - q_2*m_P

##C to P-R , scenario 5
def set_I_C_s5(e1,alfa_1,alfa_3,m_C,m_P,R,P,q_1):
    return e1*(alfa_1/m_C)*R - (alfa_3/m_P)*P - q_1



### Functions for calculating Req and Ceq in scenario 2 , Req and Peq in scenario 3 are similar
def set_R_eq_s(m_C,q_1,a_1,e1):
    return q_1*m_C/(e1*a_1)
def set_C_eq_s(m_C,r,K,q_1,a_1,e1):
    return r*(m_C/a_1)*(1- q_1*m_C/(e1*a_1*K))

#Invasibility criterion RM model

#C to R, scenario 2
def set_I_C_s2RM(e1,alfa_1,m_C,K,q_1,hC_0,q1_0):
    return alfa_1*K*(e1 - q1_0*hC_0) - q_1*m_C
#P to R, scenario 3
def set_I_P_s3RM(e_2,alfa_2,m_P,K,q_2,hP_0,q2_0):
    return alfa_2* K*(e_2 -q2_0*hP_0) - q_2*m_P

#P to C-R, scenario 4
def set_I_P_s4RM(e_2,e_3,alfa_2,alfa_3,m_P,q_2,R,C,hP_0,q2_0):
    
    return alfa_2*R*(e_2 - q2_0*hP_0) + alfa_3*C*(e_3 - q2_0*hP_0)  - q_2*m_P

#C to P-R, scenario 5
def set_I_C_s5RM(e1,e_2,alfa_1,alfa_3,m_C,m_P,R,P,q_1,t_hc,q1_0,q2_0,hP_0,hC_0):
    halfa_1 = alfa_1/m_C
    halfa_3 = alfa_3/m_P
    return (e_2/(e_2 - q2_0*hP_0)) *(halfa_1*R*(e1 - q1_0*hC_0) - q_1) - halfa_3*P*(1+t_hc*alfa_1*R)

## Functions for calculating Req and Ceq in scenario 2 , Req and Peq in scenario 3 are similar
def set_R_C_eq_sRM(m_C,r,K,q1,q1_0,alfa_1,e1,t_hc,hC_0):
    halfa_1 = alfa_1/m_C
    Req = q1 /(halfa_1*(e1 - q1_0*hC_0))
    Ceq = (1+t_hc*alfa_1*Req)/halfa_1 *r*(1-Req/K)
    return Req , Ceq

#Qualitative stability analysis using Hurwtiz determinants

##auxiliary functions
def set_D(K,a_1,a_2,a_3,e1,e_2,e_3,m_C,r):
    return K*a_1*a_2*(e1*e_3-e_2) + a_3*e_3*m_C*r

def set_b1(r,Req,K):
    return  r*Req/K

def set_b2(e1,e_2,e_3,a_1,a_2,a_3,m_C,m_P,Ceq,Req,Peq):
    return (e1*(a_1*m_P)**2*Ceq*Req + (Peq*m_C**2)*(e_2*a_2**2*Req + e_3*a_3**2*Ceq))/(m_C*m_P)**2

def set_b3(D,a_3,Ceq,Req,Peq,K,m_C,m_P):
    return a_3*Req*Ceq*Peq*D/(K*m_C*m_P**2)

##Hurwtiz criterion
def set_hdet2(b1,b2,b3):
    return b1*b2-b3


