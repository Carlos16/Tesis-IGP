from Baseclass import *
from scipy.integrate import odeint
from scipy.optimize import fsolve

##Dynamics class
class Dynamics(BSR):
    def __init__(self,workingData,initCondition,finalTime,separation,K_CP,K_RC,m_P):
        BSR.__init__(self,workingData.getParams(),workingData.getmode(),workingData.getmassLims())
        self.initCondition = initCondition
        self.finalTime = finalTime
        self.separation = separation
        self.K_CP = K_CP
        self.K_RC = K_RC
        self.m_P = m_P
        self.dR = ''
        self.dC = ''
        self.dP = ''
        self.fDict = workingData.getfDict()
        self.Kpoints = 1e4
    def setParamVals(self,K_CP,K_RC,m_P):
        """ Set the values for the three key paramaters of the model , the size ratios and the body mass """
        self.K_CP = K_CP
        self.K_RC = K_RC
        self.m_P = m_P
            
    def setinitCondition(self,initCondition):
        """ Specifices the given initial condition from which to start the simulation """
        self.initCondition = initCondition
        
    def makeinitCondition(self,case):
        
        """ se the init condition , depending on the scenario , in the Invasibility by P to C-R the initial condition
        is the equilibrium of the latter two in isolation , a similar situation is for the invasibility of C to P-R(labeled 
        scenario 5)"""
        if self.mode == "RM":
            if case == "Inv P4":
                self.initCondition[0] = self.fDict['R_eq_s2RM'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['C_eq_s2RM'](self.K_CP,self.K_RC,self.m_P)
            elif case == "Inv P5":
                self.initCondition[0] = self.fDict['R_eq_s3RM'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['P_eq_s3RM'](self.K_CP,self.K_RC,self.m_P)
        elif self.mode == "LV":
            if case == "Inv P4":
                self.initCondition[0] = self.fDict['R_eq_s2'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['C_eq_s2'](self.K_CP,self.K_RC,self.m_P)
            elif case == "Inv P5":
                self.initCondition[0] = self.fDict['R_eq_s3'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['P_eq_s3'](self.K_CP,self.K_RC,self.m_P)
                               
            
    def setfinalTime(self,finalTime):
        """ set the time until when to stop the simulation """
        self.finalTime = finalTime
    def setseparation(self,separation):
        self.separation = separation
            
    def setDynamicFunction(self):
        """ input the corresponding values for the paramters K_CP, K_RC and m_P and convert the dynamical functions dR, dC 
        and dP in three-argument functions, just depending on the value of the biomass densities R, C and P"""
             
        args = np.array([self.R,self.C,self.P,self.K_RC,self.K_CP,self.m_P])
        if self.mode == "RM":
            self.dR = lambdify((self.R,self.C,self.P),self.fDict['dRRM'](*args))
            self.dC = lambdify((self.R,self.C,self.P),self.fDict['dCRM'](*args))
            self.dP = lambdify((self.R,self.C,self.P),self.fDict['dPRM'](*args))
        else:
            self.dR = lambdify((self.R,self.C,self.P),self.fDict['dRLV'](*args))
            self.dC = lambdify((self.R,self.C,self.P),self.fDict['dCLV'](*args))
            self.dP = lambdify((self.R,self.C,self.P),self.fDict['dPLV'](*args))
            
       
    def DynamicFunction(self,X,t):
        args = np.array([X[0],X[1],X[2]])
        
        return np.array([self.dR(*args),self.dC(*args),self.dP(*args)])
    
    
        
    def Simulate(self):
        """ simulation routine using the Odeint solver from the scipy optimize package , Odeint is a python implementation
        of the ODEPACK package in FORTRAN which uses a multi-step solver in the non stiff case """
        self.setDynamicFunction()
        t = np.arange(0,self.finalTime,self.separation)
        return odeint(self.DynamicFunction,self.initCondition,t,())
    
    
    def runSimulationSimK(self,case,massLims,lowKLims,upKLims,initDirection):
        for massIndex in xrange(len(massLims)):
            
            Krange = 10**(np.linspace(lowKlims[massIndex],upKLims[massIndex],self.Kpoints))
            Kdata = []
            for K in Krange:
                self.setParamVals(K,K,massLims[massIndex])
                self.setinitCondition(case)
                Run = self.Simulate()
                Kdata.append(Run)
            WriteData(Kdata,initDirection+str(massIndex)+".csv") 


