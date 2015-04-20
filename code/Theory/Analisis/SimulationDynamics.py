from Baseclass import *
from scipy.integrate import odeint
from scipy.optimize import fsolve

##Dynamics class
class Dynamics(BSR):
    def __init__(self,workingData,initCondition,finalTime,separation,K_RC,K_CP,m_P,n=3,initMass = 1e-5):
        BSR.__init__(self,workingData.getParams(),workingData.getmode(),workingData.getxLims())
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
        self.initPop = n
        self.initMass = initMass
        self.ParamsDict = {'K_CP':self.K_CP,'K_RC':self.K_RC,'m_P':self.m_P}
        self.AssemblyInitCondition = {'Cfirst' : np.array([self.fDict['K'](self.K_RC,self.K_CP,self.m_P),self.initMass,0]) , 'Pfirst' : np.array([self.fDict['K'](self.K_RC,self.K_CP,self.m_P),0,self.initMass])}

    

    def getInitMass(self,type,n):
        if type == 'R':
            return n * self.m_P * self.K_RC * self.K_CP
        elif type == 'C':
            return n * self.m_P * self.K_CP
        else:
            return n * self.m_P


            

            


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
    
    def Bifurcation(self,focalParam,ParamRange,AssemblyType,N = 100):
        finalState =[]
        for val in ParamRange:
            self.updateFocalParam(focalParam,val)
            Run = self.AssemblySimulation(AssemblyType)[1]
            finalState.append(Run[:-N])
        return finalState

    def updateFocalParam(self,focalParam,val):
        self.ParamsDict[focalParam] = val

    def AssemblySimulation(self,type):
        # First Invasion
        initCondition = self.AssemblyInitCondition[type]
        
        self.setinitCondition(initCondition)
        
        run1 = self.Simulate()
        
        # Second Invasion
        P =Positiveformat(run1[-1])
        self.UpdateInitCondition(P,type)
        self.setinitCondition(P)
        
        run2 = self.Simulate()

        return run1,run2
        

    def UpdateInitCondition(self,P,type):
        if type == 'Cfirst':
            P[2]+=self.initMass
        else:
            P[1]+=self.initMass


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
def Positiveformat(Array):
    for i in xrange(len(Array)):
        if Array[i]< 0 :
            Array[i] = 0
            
    return Array
    
