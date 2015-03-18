# new equilibrium functions:
"""Compute equilibrium values for the 3 species taking into consideration the invasibility criterions,
i.e C will be set to 0 if inv of c to R <0 or inv of c to R-P <0 """
import numpy as np 
from numpy import linalg as LA
## R
def set_R_3(K_CP,K_RC,m_P,i_c,i_c2,eqR2,eqR3):
    return np.where(np.logical_and(i_c>0,i_c2>0),eqR2(K_CP,m_P),eqR3(K_CP,m_P))
def set_R_2(K_CP,K_RC,m_P,i_c,i_p,eqR2,eqR3,I_C2,I_P2,R_3,K):
    i_c2 = I_C2(K_CP,m_P)    
    i_p2 = I_P2(K_CP,m_P)
    return np.where(np.logical_or(np.logical_and(i_c>0,i_c2>0),np.logical_and(i_p>0,i_p2>0)),
                    R_3(K_CP,K_RC,m_P,i_c,i_c2,eqR2,eqR3),K(K_CP,m_P))

def new_R(K_CP,K_RC,m_P,f_R,I_P,I_C,I_C2,I_P2,eqR2,eqR3,R_2,R_3,K):
    i_p=I_P(K_CP,m_P)
    i_c=I_C(K_CP,m_P)
    return np.where(np.logical_and(i_p>0,i_c>0),f_R(K_CP,m_P),R_2(K_CP,K_RC,m_P,i_c,i_p,eqR2,eqR3,I_C2,I_P2,R_3,K))


##C

def set_C_2(K_CP,K_RC,m_P,i_c,eqC2,I_C2):
    i_c2 = I_C2(K_CP,m_P)
    return np.where(np.logical_and(i_c>0,i_c2>0),eqC2(K_CP,m_P),0.)
    
def new_C(K_CP,K_RC,m_P,f_C,I_P,I_C,I_C2,eqC2,C_2):
    i_p=I_P(K_CP,m_P)
    i_c=I_C(K_CP,m_P)
    return np.where(np.logical_and(i_p>0,i_c>0),f_C(K_CP,m_P),C_2(K_CP,K_RC,m_P,i_c,eqC2,I_C2))

##P
def set_P_2(K_CP,K_RC,m_P,i_p,eqP3,I_P2):
    i_p2 = I_P2(K_CP,m_P)
    return np.where(np.logical_and(i_p>0,i_p2>0),eqP3(K_CP,m_P),0.)
    
def new_P(K_CP,K_RC,m_P,f_P,I_P,I_C,I_P2,eqP3,P_2):
    i_p=I_P(K_CP,m_P)
    i_c=I_C(K_CP,m_P)
    return np.where(np.logical_and(i_p>0,i_c>0),f_P(K_CP,m_P),P_2(K_CP,K_RC,m_P,i_p,eqP3,I_P2))


def Null(*args):
    return 0
    

def setEq(args,fDict,r,mode):
    F =fDict
    if mode == 'LV':
         ReqDict={1:F['R_eq'],2:F['R_eq_s3'],3:F['R_eq_s2'],4:F['K']}
         CeqDict={1:F['C_eq'],2:Null,3:F['C_eq_s2'],4:Null} 
         PeqDict={1:F['P_eq'],2:F['P_eq_s3'],3:Null,4:Null}
        
    Req,Ceq,Peq = ReqDict[r](*args),CeqDict[r](*args),PeqDict[r](*args)
    return [Req,Ceq,Peq]
        
    
def setEigenvals(args,fDict,r,mode):
    F = fDict
    if mode == 'LV':
        EigenValsDict = {1:F['JLV'],2:F['JLVP'],3:F['JLVC']}
    
    if r<4:
         A = EigenValsDict[r](*args)
         Eigenvals = LA.eig(A)[0]
         return Eigenvals

    else:
        return [F['EigR'](*args)]

def setMTP(args,fDict,r,mode):
    F =fDict
    if mode == 'LV':
        MTPDict={1:F['MTP_C'](*args),2:2,3:2,4:1}

    return MTPDict[r]



###########Computing food chain length ##########################
def MTP_f(R_eq):
    return np.where(R_eq>0,1,0)

def MTP_O(P_eq,C_eq,R_eq,MTP_f):
    return np.where(np.logical_or(np.logical_and(C_eq>0,R_eq>0),np.logical_and(P_eq>0,R_eq>0)),2,MTP_f(R_eq))

def set_MTP(K_CP,K_RC,m_P,f_R,f_C,f_P,MTP_C,MTP_O,MTP_f,oldfR,oldfC,oldfP,I_P,I_C,I_C2,I_P2,eqR2,eqR3,eqC2,eqP3,R_2,C_2,P_2,R_3,K):
    R_eq = f_R(K_CP,K_RC,m_P,oldfR,I_P,I_C,I_C2,I_P2,eqR2,eqR3,R_2,R_3,K)
    C_eq = f_C(K_CP,K_RC,m_P,oldfC,I_P,I_C,I_C2,eqC2,C_2)
    P_eq = f_P(K_CP,K_RC,m_P,oldfP,I_P,I_C,I_P2,eqP3,P_2)
    return np.where(np.logical_and(R_eq>0,np.logical_and(C_eq>0,P_eq>0)),MTP_C(K_CP,m_P),MTP_O(P_eq,C_eq,R_eq,MTP_f))



#RM model

def setReqRM(fDict,K_CP,mP):
    Dis = fDict['Discriminant'](K_CP,mP)
    R1 = fDict['R1'](K_CP,mP)
    if Dis>0:
        bR = fDict['bR'](K_CP,mP)
        a = fDict['denR'](K_CP,mP)
        return [R1,(bR+np.sqrt(Dis))/(2*a) , (bR-np.sqrt(Dis))/(2*a)]
    elif Dis==0:
        bR = fDict['bR'](K_CP,mP)
        a =  fDict['denR'](K_CP,mP)
        return  [R1,bR/(2*a)]
    else:
        return [R1]


def setCeqRM(R,fDict,K_CP,mP):
    if R>0:
        return fDict['C_eq_RM'](R,0.,0.,K_CP,mP)
    else:
        return 0.

def setPeqRM(R,C,fDict,K_CP,mP):
    if C>0:
        return fDict['P_eq_RM'](R,0.,0.,K_CP,mP)
    else:
        return 0.


