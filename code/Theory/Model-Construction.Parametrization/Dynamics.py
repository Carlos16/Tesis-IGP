#Partial Cases

#LV 

def set_dRLVPart(R,Pred,r,K,a,m):
    return R*(r*(1-R/K) - (a/m)*Pred)

def set_dPredLV(R,Pred,a,e,m,q):
    return Pred*(e*(a/m)*R - q)

def set_dRalone(R,r,K):
    return R*r*(1-R/K)


#LV case
def set_dRLV(R,C,P,r,K,a1,a2,m_C,m_P):
    return R*(r*(1-R/K)-((a1/m_C)*C + (a2/m_P)*P))
def set_dCLV(R,C,P,a1,a3,e1,q1,m_C,m_P):
    return C*(e1*(a1/m_C)*R - (a3/m_P)*P -q1)

def set_dPLV(R,C,P,a2,a3,e2,e3,q2,m_P):
    return P*(e2*(a2/m_P)*R + e3*(a3/m_P)*C - q2)

#RM case

def set_dRRM(R,C,P,r,K,a1,a2,t_hp,t_hc,m_C,m_P):
    
    ha1 = a1/m_C
    ha2 = a2/m_P
    return R*(r*(1-R/K)-((ha1*C)/(1+t_hc*ha1*R) + (ha2*P)/(1+t_hp*ha2*R)))
def set_dCRM(R,C,P,a1,a3,e1,t_hc,t_hp,q1,m_C,m_P):
    
    ha1 = a1/m_C
    ha3 = a3/m_P

    return C*(e1*ha1*R/(1+t_hc*ha1*R) - ha3*P/(1+t_hp*ha3*C) -q1)

def set_dPRM(R,C,P,a2,a3,e2,e3,t_hp,q2,m_P):

    ha2 = a2/m_P
    ha3 = a3/m_P
    return P*(e2*ha2*R/(1+t_hp*ha2*R) + e3*ha3*C/(1+t_hp*ha3*C) - q2)

 


#Class that needs the OperationalClasses file    

