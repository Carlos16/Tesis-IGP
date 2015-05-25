
R_eq_expr='K*(-a_1*m_P*q_2 + e_3*m_C*(a_2*q_1 + a_3*r))/(K*a_1*a_2*(e_1*e_3 - e_2) + a_3*e_3*m_C*r)'
C_eq_expr='(K*a_1*a_2*e_1*m_P*q_2 - K*a_2*e_2*m_C*(a_2*q_1 + a_3*r) + a_3*m_C*m_P*q_2*r)/(a_3*(K*a_1*a_2*(e_1*e_3 - e_2) + a_3*e_3*m_C*r))'
P_eq_expr= 'm_P*(-K*a_1**2*e_1*m_P*q_2 + K*a_1*m_C*(a_2*e_2*q_1 + a_3*e_1*e_3*r) - a_3*e_3*m_C**2*q_1*r)/(a_3*m_C*(K*a_1*a_2*(e_1*e_3 - e_2) + a_3*e_3*m_C*r))'



##Equilibrium functions

def set_R_eq(K,q_1,q_2,r,a_1,a_2,a_3,e_1,e_2,e_3):
    return K*(-a_1*q_2 + e_3*(a_2*q_1 + a_3*r))/(K*a_1*a_2*(e_1*e_3 - e_2) + a_3*e_3*r)

def set_C_eq(K,q_1,q_2,r,a_1,a_2,a_3,e_1,e_2,e_3):
    return (K*a_1*a_2*e_1*q_2 - K*a_2*e_2*(a_2*q_1 + a_3*r) + a_3*q_2*r)/(a_3*(K*a_1*a_2*(e_1*e_3 - e_2) + a_3*e_3*r))

def set_P_eq(K,q_1,q_2,r,a_1,a_2,a_3,e_1,e_2,e_3):
    return (-K*a_1**2*e_1*q_2 + K*a_1*(a_2*e_2*q_1 + a_3*e_1*e_3*r) - a_3*e_3*q_1*r)/(a_3*(K*a_1*a_2*(e_1*e_3 - e_2) + a_3*e_3*r))
 
def setD(K,a_1,a_2,a_3,e_1,e_2,e_3,r):
    return K*a_1*a_2*(e_1*e_3  - e_2) + a_3*e_3*r

def setDBound(K,a_1,a_2,a_3,e_1,e_2,e_3,m_C,r):
    return a_3*e_3*r/(K*a_1*a_2)



def setEqPNum_RM(K,q_1,q_2,r,a1,a2,a3,e_1,e_2,e_3,thc,thp,R,C,P,mP,mC,q2_0,q1_0,hC_0,hP_0):
    return -mP*(R**2*a1*a2*a3*e_1*e_2*thp - R**2*a1*a2*a3*e_1*hP_0*q2_0*thp - R**2*a1*a2*a3*e_2*mC*q_1*thc*thp + R**2*a1*a2*a3*hP_0*mC*q2_0*q_1*thc*thp - R**2*a1*a2*e_1*e_3*thp + R**2*a1*a2*e_1*hP_0*q2_0*thp + R**2*a1*a2*e_3*mC*q_1*thc*thp - R**2*a1*a2*hP_0*mC*q2_0*q_1*thc*thp - R*a1*a3*e_1*mP*q_2*thp + R*a1*a3*mC*mP*q_1*q_2*thc*thp - R*a1*e_1*e_3 + R*a1*e_1*hP_0*q2_0 + R*a1*e_3*mC*q_1*thc - R*a1*hP_0*mC*q2_0*q_1*thc - R*a2*a3*e_2*mC*q_1*thp + R*a2*a3*hP_0*mC*q2_0*q_1*thp + R*a2*e_3*mC*q_1*thp - R*a2*hP_0*mC*q2_0*q_1*thp + a3*mC*mP*q_1*q_2*thp + e_3*mC*q_1 - hP_0*mC*q2_0*q_1)

def setEqPDen_RM(K,q_1,q_2,r,a1,a2,a3,e_1,e_2,e_3,thc,thp,R,C,P,mP,mC,q2_0,q1_0,hC_0,hP_0):
    return a3*mC*(R*a1*e_3*thc - R*a1*hP_0*q2_0*thc + e_3 - hP_0*q2_0)



def setEqCNum_RM(q_2,mP,a2,R,e_2,q2_0,hP_0):
    return q_2*mP - a2*R*(e_2 - q2_0*hP_0)

def setEqCDen_RM(e_3,q2_0,hP_0):
    return e_3 - q2_0*hP_0

def setDis(K,q_1,q_2,r,a1,a2,a3,e_1,e_2,e_3,thc,thp,mP,mC,q2_0,q1_0,hC_0,hP_0):
    return K**2*a1**2*a2**2*a3**2*e_2**2 - 2*K**2*a1**2*a2**2*a3**2*e_2*hP_0*q2_0 + K**2*a1**2*a2**2*a3**2*hP_0**2*q2_0**2 - 2*K**2*a1**2*a2**2*a3*e_1*e_2*e_3 + 2*K**2*a1**2*a2**2*a3*e_1*e_2*hP_0*q2_0 + 2*K**2*a1**2*a2**2*a3*e_1*e_3*hP_0*q2_0 - 2*K**2*a1**2*a2**2*a3*e_1*hP_0**2*q2_0**2 + 2*K**2*a1**2*a2**2*a3*e_2*e_3*mC*q_1*thc - 2*K**2*a1**2*a2**2*a3*e_2*hP_0*mC*q2_0*q_1*thc - 2*K**2*a1**2*a2**2*a3*e_3*hP_0*mC*q2_0*q_1*thc + 2*K**2*a1**2*a2**2*a3*hP_0**2*mC*q2_0**2*q_1*thc + K**2*a1**2*a2**2*e_1**2*e_3**2 - 2*K**2*a1**2*a2**2*e_1**2*e_3*hP_0*q2_0 + K**2*a1**2*a2**2*e_1**2*hP_0**2*q2_0**2 - 2*K**2*a1**2*a2**2*e_1*e_3**2*mC*q_1*thc + 4*K**2*a1**2*a2**2*e_1*e_3*hP_0*mC*q2_0*q_1*thc - 2*K**2*a1**2*a2**2*e_1*hP_0**2*mC*q2_0**2*q_1*thc + K**2*a1**2*a2**2*e_3**2*mC**2*q_1**2*thc**2 - 2*K**2*a1**2*a2**2*e_3*hP_0*mC**2*q2_0*q_1**2*thc**2 + K**2*a1**2*a2**2*hP_0**2*mC**2*q2_0**2*q_1**2*thc**2 + 2*K**2*a1**2*a2*a3**2*e_2*e_3*mC*r*thc - 2*K**2*a1**2*a2*a3**2*e_2*hP_0*mC*q2_0*r*thc - 2*K**2*a1**2*a2*a3**2*e_3*hP_0*mC*q2_0*r*thc + 2*K**2*a1**2*a2*a3**2*hP_0**2*mC*q2_0**2*r*thc - 2*K**2*a1**2*a2*a3*e_1*e_3**2*mC*r*thc + 4*K**2*a1**2*a2*a3*e_1*e_3*hP_0*mC*q2_0*r*thc - 2*K**2*a1**2*a2*a3*e_1*hP_0**2*mC*q2_0**2*r*thc + 2*K**2*a1**2*a2*a3*e_3**2*mC**2*q_1*r*thc**2 - 4*K**2*a1**2*a2*a3*e_3*hP_0*mC**2*q2_0*q_1*r*thc**2 + 2*K**2*a1**2*a2*a3*hP_0**2*mC**2*q2_0**2*q_1*r*thc**2 + K**2*a1**2*a3**2*e_3**2*mC**2*r**2*thc**2 - 2*K**2*a1**2*a3**2*e_3*hP_0*mC**2*q2_0*r**2*thc**2 + K**2*a1**2*a3**2*hP_0**2*mC**2*q2_0**2*r**2*thc**2 - 4*K*a1**2*a3**2*e_3*mC*mP*q_2*r*thc + 4*K*a1**2*a3**2*hP_0*mC*mP*q2_0*q_2*r*thc - 2*K*a1*a2*a3**2*e_2*e_3*mC*r + 2*K*a1*a2*a3**2*e_2*hP_0*mC*q2_0*r + 2*K*a1*a2*a3**2*e_3*hP_0*mC*q2_0*r - 2*K*a1*a2*a3**2*hP_0**2*mC*q2_0**2*r + 2*K*a1*a2*a3*e_1*e_3**2*mC*r - 4*K*a1*a2*a3*e_1*e_3*hP_0*mC*q2_0*r + 2*K*a1*a2*a3*e_1*hP_0**2*mC*q2_0**2*r + 2*K*a1*a2*a3*e_3**2*mC**2*q_1*r*thc - 4*K*a1*a2*a3*e_3*hP_0*mC**2*q2_0*q_1*r*thc + 2*K*a1*a2*a3*hP_0**2*mC**2*q2_0**2*q_1*r*thc + 2*K*a1*a3**2*e_3**2*mC**2*r**2*thc - 4*K*a1*a3**2*e_3*hP_0*mC**2*q2_0*r**2*thc + 2*K*a1*a3**2*hP_0**2*mC**2*q2_0**2*r**2*thc + a3**2*e_3**2*mC**2*r**2 - 2*a3**2*e_3*hP_0*mC**2*q2_0*r**2 + a3**2*hP_0**2*mC**2*q2_0**2*r**2

def setden_R(K,q_1,q_2,r,a1,a2,a3,e_1,e_2,e_3,thc,thp,mP,mC,q2_0,q1_0,hC_0,hP_0):
    return 2*a1*a3*mC*r*thc*(e_3 - hP_0*q2_0)

def setb_R(K,q_1,q_2,r,a1,a2,a3,e_1,e_2,e_3,thc,thp,mP,mC,q2_0,q1_0,hC_0,hP_0):
    return K*a1*a2*a3*e_2 - K*a1*a2*a3*hP_0*q2_0 - K*a1*a2*e_1*e_3 + K*a1*a2*e_1*hP_0*q2_0 + K*a1*a2*e_3*mC*q_1*thc - K*a1*a2*hP_0*mC*q2_0*q_1*thc + K*a1*a3*e_3*mC*r*thc - K*a1*a3*hP_0*mC*q2_0*r*thc - a3*e_3*mC*r + a3*hP_0*mC*q2_0*r


def setRoot1(K,q_1,q_2,r,a1,a2,a3,e_1,e_2,e_3,thc,thp,mP,mC,q2_0,q1_0,hC_0,hP_0):
    return (a3*hP_0*q2_0 + e_3 - hP_0*q2_0)/(a2*thp*(a3*e_2 - a3*hP_0*q2_0 - e_3 + hP_0*q2_0))





