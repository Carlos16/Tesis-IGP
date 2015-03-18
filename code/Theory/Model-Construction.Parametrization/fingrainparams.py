from coarsegrainparams import *
from inva_fcl_stab import *
from Eq import *
from Dynamics import *
from sympy import Matrix,sqrt

def construct_param_dict(params,K_RC,K_CP,m_P):
    ###scaling constants
    w=params['w']
    pd=params['pd'] # in 3D and 0.21 in 2D
    pv=params['pv']
    Er=params['Er'] ;Ek=params['Ek']
    ER=params['ER'];EC=params['EC'];EP=params['EP'];
    Eq1=params['Eq1'];Eq2=params['Eq2']


    #capture success function
    a = params['a']
    b = params['b']
    c = params['c']
    formC = params['formC']
    formPC = params['formPC']
    formPR = params['formPR']
    
    ###variables
    TR= params['TR'] ;TC= params['TC'];TP=params['TP'];D_R= params['D_R']; D_C= params['D_C']
    K_RP=K_RC*K_CP
    fmC=params['fmC'];thermyR=params['thermyR']
    thermyC=params['thermyC'];thermyP=params['thermyP']
    fmPR=params['fmPR']
    fmPC=params['fmPC']
    m_C = K_CP*m_P;m_R = K_RP*m_P
    ###normalization constants and boltzmann constant
    r0 = params['r0']
    k0 = params['k0'] # will depend on the productivity of the habitat
    alfa01 = alfa02 = params['alfa012']  # will depedend on the dimension of the habitat 
    alfa03 = params['alfa03']
    d0= params['d0']
    q10 = params['q10'];q20 = params['q20'];
    v0R = params['v0R'];v0C =params['v0C'];v0P =params['v0P'];k = b_k
    t_h0C = params['t_h0C'];t_h0P = params['t_h0P'] 

    ##efficience values, later check 
    ##for the sensibility of the results to changes in it
    e_1= params['e_1'];e_2=params['e_2'];e_3=params['e_3']
    
    #intrapopulation parameters
    q1=set_q1(q10,m_C,w,Eq1,TR,k)
    q2=set_q2(q20,m_P,w,Eq2,TC,k)
    K=set_K(k0,m_R,w,Ek,TR,k)
    r=set_r(r0,m_R,w,Er,TR,k)

    #interpopulation parameters
    alfa1=set_alfa(m_C,alfa01,K_RC,pv,pd,TR,TC,ER,EC,D_R,v0R,v0C,g,alfa,fmC,thermyR,thermyC,k,a,b,c,formC)
    alfa2=set_alfa(m_P,alfa02,K_RP,pv,pd,TR,TP,ER,EP,D_R,v0R,v0P,g,alfa,fmPR,thermyR,thermyP,k,a,b,c,formPR)
    alfa3=set_alfa(m_P,alfa03,K_CP,pv,pd,TC,TP,EC,EP,D_C,v0C,v0P,g,alfa,fmPC,thermyC,thermyP,k,a,b,c,formPC)

    t_hp = set_th(t_h0P,m_P,w,EP,k,TP)
    t_hc = set_th(t_h0C,m_C,w,EC,k,TC)
    param_dict={'q1':q1,'q2':q2,'K':K,'r':r,'alfa1':alfa1,'alfa2':alfa2,'alfa3':alfa3,'t_hp':t_hp,'t_hc':t_hc}
       
    return param_dict

def construct_equilibrium(params,par_dict,K_RC,K_CP,m_P):

    #intrapopulation parameters
    q1=par_dict['q1']
    q2=par_dict['q2']
    q1_0 = params['q10']
    q2_0 = params['q20']
    hC_0 = params['t_h0C']
    hP_0 = params['t_h0P']
    
    K=par_dict['K']
    r=par_dict['r']
    
    m_C = K_CP*m_P

    #interpopulation parameters
    alfa1=par_dict['alfa1']
    alfa2=par_dict['alfa2']
    alfa3=par_dict['alfa3']
    t_hc = par_dict['t_hc']
    t_hp = par_dict['t_hp']
    e_1=params['e_1']
    e_2=params['e_2']
    e_3=params['e_3']
    

    ## Equilibrium values
    ###Sc2
    R_eq_s2 = set_R_eq_s(m_C,q1,alfa1,e_1)
    C_eq_s2 = set_C_eq_s(m_C,r,K,q1,alfa1,e_1)

    R_eq_s2RM, C_eq_s2RM = set_R_C_eq_sRM(m_C,r,K,q1,q1_0,alfa1,e_1,t_hc,hC_0)

    ###Sc3
    R_eq_s3 = set_R_eq_s(m_P,q2,alfa2,e_2)
    P_eq_s3 = set_C_eq_s(m_P,r,K,q2,alfa2,e_2)

    R_eq_s3RM , P_eq_s3RM = set_R_C_eq_sRM(m_P,r,K,q2,q2_0,alfa2,e_2,t_hp,hP_0)
    ###full system
    R_eq = set_R_eq(m_P,m_C,K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3)
    C_eq = set_C_eq(m_P,m_C,K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3)
    P_eq = set_P_eq(m_P,m_C,K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3)
    
    D = setD(K,alfa1,alfa2,alfa3,e_1,e_2,e_3,m_C,r)
    DBound= setDBound(K,alfa1,alfa2,alfa3,e_1,e_2,e_3,m_C,r)
    #Roots for Req
    R1 = setRoot1(K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3,t_hc,t_hp,m_P,m_C,q2_0,q1_0,hC_0,hP_0)
    Dis = setDis(K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3,t_hc,t_hp,m_P,m_C,q2_0,q1_0,hC_0,hP_0)
    bR = setb_R(K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3,t_hc,t_hp,m_P,m_C,q2_0,q1_0,hC_0,hP_0)
    denR = setden_R(K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3,t_hc,t_hp,m_P,m_C,q2_0,q1_0,hC_0,hP_0)

    R2 = (bR + sqrt(Dis))/(2*denR)
    R3 = (bR - sqrt(Dis))/(2*denR)
    
    
    
    eq_dict={'R_eq_s2':R_eq_s2,'C_eq_s2':C_eq_s2,'R_eq_s3':R_eq_s3,'P_eq_s3':P_eq_s3,'R_eq':R_eq,'C_eq':C_eq,'P_eq':P_eq,
             'R_eq_s2RM':R_eq_s2RM,'C_eq_s2RM':C_eq_s2RM,'R_eq_s3RM':R_eq_s3RM,'P_eq_s3RM':P_eq_s3RM,'R1':R1,'Discriminant':Dis,'R2':R2,'R3':R3,'bR':bR,'denR':denR,'D' : D,'DBound':DBound}
    return eq_dict



def construct_inv_boundaries(params,par_dict,eq_dict,K_RC,K_CP,m_P):
    #intrapop params
    q1=par_dict['q1']
    q2=par_dict['q2']
    K =par_dict['K']
    m_C= K_CP*m_P
    q1_0 = params['q10']
    q2_0 = params['q20']
    hC_0 = params['t_h0C']
    hP_0 = params['t_h0P']

    #interpop params
    alfa1=par_dict['alfa1']
    alfa2=par_dict['alfa2']
    alfa3=par_dict['alfa3']
    e_1=params['e_1']
    e_2=params['e_2']
    e_3=params['e_3']
    

    t_hc = par_dict['t_hc']
    t_hp = par_dict['t_hp']


    #eq values

    #L-V
    R_eq_s2 = eq_dict['R_eq_s2']
    C_eq_s2 = eq_dict['C_eq_s2']
    P_eq_s3 = eq_dict['P_eq_s3']
    R_eq_s3 = eq_dict['R_eq_s3']
    #R-M
    R_eq_s2RM = eq_dict['R_eq_s2RM']
    C_eq_s2RM = eq_dict['C_eq_s2RM']
    R_eq_s3RM = eq_dict['R_eq_s3RM']
    P_eq_s3RM = eq_dict['P_eq_s3RM']
    
    ##Invasibility boundaries

    #L-V
    I_C_s2 = set_I_C_s2(e_1,alfa1,m_C,K,q1)
    I_P_s3 = set_I_P_s3(e_2,alfa2,m_P,K,q2)
    I_P_s4 = set_I_P_s4(e_2,e_3,alfa2,alfa3,m_P,q2,R_eq_s2,C_eq_s2)
    I_C_s5 = set_I_C_s5(e_1,alfa1,alfa3,m_C,m_P,R_eq_s3,P_eq_s3,q1)


    
    #R-M
    I_C_s2RM = set_I_C_s2RM(e_1,alfa1,m_C,K,q1,hC_0,q1_0)
    I_P_s3RM = set_I_P_s3RM(e_2,alfa2,m_P,K,q2,hP_0,q2_0)
    I_P_s4RM = set_I_P_s4RM(e_2,e_3,alfa2,alfa3,m_P,q2,R_eq_s2RM,C_eq_s2RM,hP_0,q2_0)
    I_C_s5RM = set_I_C_s5RM(e_1,e_2,alfa1,alfa3,m_C,m_P,R_eq_s3RM,P_eq_s3RM,q1,t_hc,q1_0,q2_0,hP_0,hC_0)    

    inv_dict= {'I_C_s2':I_C_s2,'I_P_s3':I_P_s3,'I_P_s4':I_P_s4,'I_C_s5':I_C_s5,
               'I_C_s2RM':I_C_s2RM,'I_P_s3RM':I_P_s3RM,'I_P_s4RM':I_P_s4RM,'I_C_s5RM':I_C_s5RM}

    return inv_dict

def Trophic_position(params,par_dict,eq_dict,m_P):
    
    R_eq = eq_dict['R_eq']
    alfa2 = par_dict['alfa2']
    q2 = par_dict['q2']
    e_2 = params['e_2']
    
    #Trophic position in the coexistence domain
    MTP_C= set_MTP_C(R_eq,alfa2,m_P,q2,e_2)
    return MTP_C
    
def Stability(params,par_dict,eq_dict,K_RC,K_CP,m_P):
    #intrapop params
    K=par_dict['K']
    r=par_dict['r']
    m_C = K_CP*m_P
    #interpop params
    alfa1=par_dict['alfa1']
    alfa2=par_dict['alfa2']
    alfa3=par_dict['alfa3']
    e_1=params['e_1']
    e_2=params['e_2']
    e_3=params['e_3']
    #equilibrium
    R_eq= eq_dict['R_eq']
    C_eq = eq_dict['C_eq']
    P_eq = eq_dict['P_eq']
    
    
    ##Stability
    D = set_D(K,alfa1,alfa2,alfa3,e_1,e_2,e_3,m_C,r)
    a1 = set_a1(r,R_eq,K)
    a2 = set_a2(e_1,e_2,e_3,alfa1,alfa2,alfa3,m_C,m_P,C_eq,R_eq,P_eq)
    a3 = set_a3(D,alfa3,C_eq,R_eq,P_eq,K,m_C,m_P)
    hd2 = set_hdet2(a1,a2,a3)

    return hd2

def Jacobian(dR,dC,dP,R,C,P):
    X = Matrix([dR,dC,dP])
    Y = Matrix([R,C,P])
    return X.jacobian(Y)

def Jacobian2(dX,dY,X,Y):
    A = Matrix([dX,dY])
    B = Matrix([X,Y])
    return A.jacobian(B)

def setJacobianDict(DynamicsDict,R,C,P):
    dRLV = DynamicsDict['dRLV']
    dCLV = DynamicsDict['dCLV']
    dPLV = DynamicsDict['dPLV']
    dRRM = DynamicsDict['dRRM']
    dCRM = DynamicsDict['dCRM']
    dPRM = DynamicsDict['dPRM']
    
    dRLVP = DynamicsDict['dRLVP']
    dRLVC = DynamicsDict['dRLVC']
    dPLVP = DynamicsDict['dPLVP']
    dCLVC = DynamicsDict['dCLVC']
   
    JLV  = Jacobian(dRLV,dCLV,dPLV,R,C,P)
    JRM = Jacobian(dRRM,dCRM,dPRM,R,C,P)
    JLVP = Jacobian2(dRLVP,dPLVP,R,P)
    JLVC = Jacobian2(dRLVC,dCLVC,R,C)
    
    return {'JLV':JLV,'JRM':JRM,'JLVP':JLVP,'JLVC':JLVC}

def ConstructDynamicalFunctions(params,par_dict,K_RC,K_CP,m_P,R,C,P):
    #intrapopulation parameters
    q1=par_dict['q1']
    q2=par_dict['q2']
    K=par_dict['K']
    r=par_dict['r']
    e_1=params['e_1']
    e_2=params['e_2']
    e_3=params['e_3']
    m_C = K_CP*m_P
    q2_0 = params['q20']
    q1_0 = params['q10']
    #interpopulation parameters
    alfa1=par_dict['alfa1']
    alfa2=par_dict['alfa2']
    alfa3=par_dict['alfa3']
    t_hp=par_dict['t_hp']
    t_hc=par_dict['t_hc']
    hC_0=params['t_h0C']
    hP_0=params['t_h0P']
    
    #Construct LV functions
    dRLV=set_dRLV(R,C,P,r,K,alfa1,alfa2,m_C,m_P)
    dPLV=set_dPLV(R,C,P,alfa2,alfa3,e_2,e_3,q2,m_P)
    dCLV=set_dCLV(R,C,P,alfa1,alfa3,e_1,q1,m_C,m_P)

    dRLVP = set_dRLVPart(R,P,r,K,alfa2,m_P)
    dPLVP = set_dPredLV(R,P,alfa2,e_2,m_P,q2)

    dRLVC = set_dRLVPart(R,C,r,K,alfa1,m_C)
    dCLVC = set_dPredLV(R,C,alfa1,e_1,m_C,q1)

        
    #Construct RM functions
    dRRM = set_dRRM(R,C,P,r,K,alfa1,alfa2,t_hp,t_hc,m_C,m_P)
    dPRM = set_dPRM(R,C,P,alfa2,alfa3,e_2,e_3,t_hp,q2,m_P) 
    dCRM = set_dCRM(R,C,P,alfa1,alfa3,e_1,t_hc,t_hp,q1,m_C,m_P)


    #RM eq expresions

    CNum_eq_RM = setEqCNum_RM(q2,m_P,alfa2,R,e_2,q2_0,hP_0)
    CDen_eq_RM = setEqCDen_RM(e_3,q2_0,hP_0)
    PNum_eq_RM = setEqPNum_RM(K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3,t_hc,t_hp,R,C,P,m_P,m_C,q2_0,q1_0,hC_0,hP_0)
    PDen_eq_RM = setEqPDen_RM(K,q1,q2,r,alfa1,alfa2,alfa3,e_1,e_2,e_3,t_hc,t_hp,R,C,P,m_P,m_C,q2_0,q1_0,hC_0,hP_0)


    C_eq_RM = CNum_eq_RM/CDen_eq_RM
    P_eq_RM = PNum_eq_RM/PDen_eq_RM

    DynamicsDict={'dRLV':dRLV,'dPLV':dPLV,'dCLV':dCLV,'dRRM':dRRM,'dPRM':dPRM,'dCRM':dCRM,'C_eq_RM':C_eq_RM,'P_eq_RM':P_eq_RM,'PNum_eq_RM':PNum_eq_RM,'CNum_eq_RM':CNum_eq_RM,'dRLVP':dRLVP,'dPLVP':dPLVP,'dRLVC':dRLVC,'dCLVC':dCLVC,'EigR':-r}
    return DynamicsDict

    


    

    
    
