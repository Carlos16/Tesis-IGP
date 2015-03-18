#Computing food chain length
def MTP_f(R_eq):
    return np.where(R_eq>0,1,0)

def MTP_O(P_eq,C_eq,R_eq,MTP_f):
    return np.where(np.logical_or(np.logical_and(C_eq>0,R_eq>0),np.logical_and(P_eq>0,R_eq>0)),2,MTP_f(R_eq))

def set_MTP(K_CP,K_RC,m_P,f_R,f_C,f_P,MTP_C,MTP_O,MTP_f):
    R_eq = f_R(K_CP,K_RC,m_P)
    C_eq = f_C(K_CP,K_RC,m_P)
    P_eq = f_P(K_CP,K_RC,m_P)
    return np.where(np.logical_and(R_eq>0,np.logical_and(C_eq>0,P_eq>0)),MTP_C(K_CP,K_RC,m_P),MTP_O(P_eq,C_eq,R_eq,MTP_f))
