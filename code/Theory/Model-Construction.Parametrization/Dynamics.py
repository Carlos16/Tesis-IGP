#Partial Cases

##LV 

def set_dRLVPart(R,Pred,r,K,a):
    return R*(r*(1-R/K) - a*Pred)

def set_dPredLV(R,Pred,a1,e,q):
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

def set_dRLV(R,C,P,r,K,a1,a2):
    return R*(r*(1-R/K)-(a1*C + a2 *P))
def set_dCLV(R,C,P,a1,a3,e1,q1):
    return C*(e1*a1*R - a3 *P -q1)

def set_dPLV(R,C,P,a2,a3,e2,e3,q2):
    return P*(e2*a2*R + e3*a3*C - q2)


##LV case adimensional

def set_LVAdim(x,y,z,r,K,a1,a2,a3,e1,e2,e3,q1,q2):
    w1 = a1*K/r
    w2 = a2*K/r
    w3 = a3*K/r
    d1 = q1/r
    d2 = q2/r
    
    dx = x*((1-x)- (w1*y+ w2*z))
    dy = y*(e1*w1*x - w3*z - d1)
    dz = z*(e2*w2*x + e3*w3*y - d2)

    return dx,dy,dz


           
    
def set_IsoclinesLVAdim(x,y,z,r,K,a1,a2,a3,e1,e2,e3,q1,q2):
    w1 = a1*K/r
    w2 = a2*K/r
    w3 = a3*K/r
    d1 = q1/r
    d2 = q2/r

    R,C,P = IsoclinesLV(y,z,w1,w2,w3,e1,e2,e3,d1,d2)

    return R,C,P
    
def IsoclinesLV(y,z,w1,w2,w3,e1,e2,e3,d1,d2):
    R = 1 - w1*y - w2*z
    C = (w3* z  + d1)/(e1*w1)
    P = (-e3*w3*y + d2)/(e2*w2)
    
    return R,C,P


    

##RM case

def set_dRRM(R,C,P,r,K,a1,a2,a3,t_hp,t_hc,m_C,m_P):
    
    return R*(r*(1-R/K)-((a1*C)/(1+t_hc*a1*m_C*R) + (a2*P)/(1+t_hp*m_P*(a2*R+a3*C))))
def set_dCRM(R,C,P,a1,a2,a3,e1,t_hc,t_hp,q1,m_C,m_P):
    
    return C*(e1*a1*R/(1+t_hc*m_C*a1*R) - a3*P/(1+t_hp*m_P*(a3*C+a2*R)) -q1)

def set_dPRM(R,C,P,a2,a3,e2,e3,t_hp,q2,m_P):

    return P*((e2*a2*R+e3*a3*C)/(1+t_hp*m_P*(a3*C+a2*R)) - q2)

 
