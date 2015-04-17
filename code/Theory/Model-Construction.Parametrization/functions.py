from fingrainparams import *
from sympy import symbols,lambdify

def func_transform(params,K_CP,K_RC,m_P,R,C,P,sim = True,bottom = False):
    
    if sim :
        args = (K_CP,K_CP,m_P)
    else:
        args = (K_RC,K_CP,m_P)
    
    # Construct dictionaries of expressions
    args1 = (params,)+args
    
    par_dict=construct_param_dict(*args1)
    
    args2 = (params,par_dict)+args

    eq_dict = construct_equilibrium(*args2)

    args3 = (params,par_dict,eq_dict)+args

    inv_dict= construct_inv_boundaries(*args3)
    MTP_C= Trophic_position(params,par_dict,eq_dict,m_P)
    hd2= Stability(*args3)

    args4 = args2 + (R,C,P)

    DynamicsDict = ConstructDynamicalFunctions(*args4)
    JacobianDict = setJacobianDict(DynamicsDict,R,C,P)
    
    #Initial Setting#
    f_dict_keys = ['hd2','MTP_C']
    f_dict_vals = [hd2,MTP_C]

    #Load Keys#
    if bottom:
        m_R = symbols('m_R')
        f_dict = { f_dict_keys[i] : lambdify(setArgs(args),f_dict_vals[i]) for i in range(len(f_dict_keys))}
        inv_dict = {k:lambdify(setArgs(args),v) for (k,v) in inv_dict.items()}
        eq_dict = {k:lambdify(setArgs(args),v) for (k,v) in eq_dict.items()}
        DynamicsDict = {k:lambdify((R,C,P)+tuple(setArgs(args)),v) for (k,v) in DynamicsDict.items()}
        JacobianDict = {k:lambdify((R,C,P,K_CP,K_RC,m_R),v) for (k,v) in JacobianDict.items()}
    else:
        f_dict= { f_dict_keys[i] : lambdify(setArgs(args),f_dict_vals[i]) for i in range(len(f_dict_keys))}
        param_dict={k:lambdify(setArgs(args),v) for (k,v) in par_dict.items()}
        inv_dict = {k:lambdify(setArgs(args),v) for (k,v) in inv_dict.items()}
        eq_dict = {k:lambdify(setArgs(args),v) for (k,v) in eq_dict.items()}
        DynamicsDict = {k:lambdify((R,C,P)+tuple(setArgs(args)),v) for (k,v) in DynamicsDict.items()}
        JacobianDict = {k:lambdify((R,C,P)+tuple(setArgs(args)),v) for (k,v) in JacobianDict.items()}
        

    f_dict = dict(f_dict.items() + inv_dict.items() + param_dict.items()+eq_dict.items()+ DynamicsDict.items() + JacobianDict.items())
    

    return f_dict

def setArgs(args):
    n1 = len(args)
    n2 = len(set(args))
    if n2<n1:
        return (args[0],args[-1])
    else:
        return args



    
