#Partial Cases

##LV 

def set_dRLVPart(R,Pred,r,K,a,m):
    a1 = a/m
    return R*(r*(1-R/K) - a1*Pred)

def set_dPredLV(R,Pred,a,e,m,q):
    a1 = a/m
    return Pred*(e*a1*R - q)

def set_dRalone(R,r,K):
    return R*r*(1-R/K)

##RM


def set_dRRMPart(R,Pred,r,K,a,m,t_h):
    a1 = a/m
    return R*(r*(1-R/K) - a1*Pred/(1+t_h*a*R))

def set_dPredRM(R,Pred,r,K,a,e,m,q,t_h):
    a1 = a/m
    return Pred*(e*a1*R/(1+t_h*a*R) - q)

#Complete System
##LV case

def set_dRLV(R,C,P,r,K,a1,a2,m_C,m_P):
    return R*(r*(1-R/K)-((a1/m_C)*C + (a2/m_P)*P))
def set_dCLV(R,C,P,a1,a3,e1,q1,m_C,m_P):
    return C*(e1*(a1/m_C)*R - (a3/m_P)*P -q1)

def set_dPLV(R,C,P,a2,a3,e2,e3,q2,m_P):
    return P*(e2*(a2/m_P)*R + e3*(a3/m_P)*C - q2)

##RM case

def set_dRRM(R,C,P,r,K,a1,a2,a3,t_hp,t_hc,m_C,m_P):
    
    return R*(r*(1-R/K)-((a1*C)/(1+t_hc*a1*m_C*R) + (a2*P)/(1+t_hp*m_P*(a2*R+a3*C))))
def set_dCRM(R,C,P,a1,a2,a3,e1,t_hc,t_hp,q1,m_C,m_P):
    
    return C*(e1*a1*R/(1+t_hc*m_C*a1*R) - a3*P/(1+t_hp*m_P*(a3*C+a2*R)) -q1)

def set_dPRM(R,C,P,a2,a3,e2,e3,t_hp,q2,m_P):

    return P*((e2*a2*R+e3*a3*C)/(1+t_hp*m_P*(a3*C+a2*R)) - q2)

 
