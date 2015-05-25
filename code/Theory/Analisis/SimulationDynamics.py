from Baseclass import *
from scipy.integrate import odeint
from scipy.integrate import ode
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
        self.Jacobian = ''
        self.fDict = workingData.getfDict()
        self.Kpoints = 1e4
        self.initPop = n
        self.initMass = initMass
        self.ParamsDict = {'K_CP':self.K_CP,'K_RC':self.K_RC,'m_P':self.m_P}
        self.AssemblyInitCondition = {'Cfirst' : np.array([1.,self.initMass,0]) , 'Pfirst' : np.array([1.,0,self.initMass])}

    def getIsoclines(self,Y,Z):
        R = self.fDict['RIsoLVa'](1.,Y,Z,self.K_RC,self.K_CP,self.m_P)
        C = self.fDict['CIsoLVa'](1.,Y,Z,self.K_RC,self.K_CP,self.m_P)
        P = self.fDict['PIsoLVa'](1.,Y,Z,self.K_RC,self.K_CP,self.m_P)

        return R,C,P
        

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
            self.Jacobian = lambdify((self.R,self.C,self.P),self.fDict['JRM'](*args))
        else:
            dx= self.fDict['dxLVa'](*args)
            dy= self.fDict['dyLVa'](*args)
            dz= self.fDict['dzLVa'](*args)

            self.dR = lambdify((self.R,self.C,self.P),dx)
            self.dC = lambdify((self.R,self.C,self.P),dy)
            self.dP = lambdify((self.R,self.C,self.P),dz)           

            self.Jacobian = lambdify((self.R,self.C,self.P),Jacobian(dx,dy,dz,self.R,self.C,self.P))
            
       
    def DynamicFunction2(self,t,Y):
        return np.array([self.dR(*Y),self.dC(*Y),self.dP(*Y)])

    def DynamicFunction(self,X,t):
        
        return np.array([self.dR(*X),self.dC(*X),self.dP(*X)])
        
    def JacobianFunction(self,t,X):
        return np.array(self.Jacobian(*X))

    def JacobianFunction2(self,X,t):
        return np.array(self.Jacobian(*X))
        
        
    def Simulate(self):
        """ simulation routine using the Odeint solver from the scipy optimize package , Odeint is a python implementation
        of the ODEPACK package in FORTRAN which uses a multi-step solver in the non stiff case """
        self.setDynamicFunction()
        t = np.arange(0,self.finalTime,self.separation)
        return odeint(self.DynamicFunction,self.initCondition,t,(),Dfun= self.JacobianFunction2,full_output=1)

    def Simulate2(self,integrator='lsoda', t0=0,N=3):
        self.setDynamicFunction()
        tf = self.finalTime
        dt = self.separation
        r = ode(self.DynamicFunction2,jac=self.JacobianFunction).set_integrator(integrator)

        r.set_initial_value(self.initCondition,t0)

        FinalS=[]

        while r.successful() and r.t < (tf-dt*N):
            r.integrate(r.t+dt)

        while r.t < tf:
            r.integrate(r.t+dt)
            if r.successful():
                FinalS.append(r.y)
            else:
                break
                
        return FinalS,r        
        
        
    
    def BifurcationLV(self,focalParam,ParamRange,AssemblyType,N = 10):
        finalState =[]
        for val in ParamRange:
            self.updateFocalParam(focalParam,val)
            I1,I2,I3,I4 = self.getInv()
            if AssemblyType == 'Cfirst':
                if I1 > 0 and I3 <0 :
                    X = self.getEqScenario(2)                    
                    finalState.append([X])
                elif I1<0 and I2>0:
                    X = self.getEqScenario(3)
                    finalState.append([X])
                elif I1<0 and I2<0:
                    finalState.append([np.array([1.,0.,0.])])
                else:
                    Run = self.AssemblySimulation(AssemblyType)
                    finalState.append(Run)
                
            if AssemblyType == 'Pfirst':

                if I2>0 and I4 <0 :
                    X = self.getEqScenario(3)
                    finalState.append([np.array([xf,yf,zf])])
                elif I2 <0 and I1>0:
                    X = self.getEqScenario(2)
                    finalState.append([X])
                elif I1<0 and I2<0:
                    finalsState.append([np.array([1.,0.,0.])])
                else:
                    Run = self.AssemblySimulation(AssemblyType)
                    finalState.append(Run)
        return finalState

    def setAndWriteBifurcation(self,focalParam,ParamRange,AssemblyType,Direction,N = 10, header = ['x','y','z']):
        F = self.BifurcationLV(focalParam,ParamRange,AssemblyType,N)
        footer = ConstructFooter(self.params)
        SaveData(F,direction,header,footer)
        
    def getEqScenario(self,n):

        K = self.fDict['K'](self.K_RC,self.K_CP,self.m_P)
        if n == 2:
            x = self.fDict['R_eq_s2'](self.K_RC,self.K_CP,self.m_P) /K
            y = self.fDict['C_eq_s2'](self.K_RC,self.K_CP,self.m_P) /K
            z = 0.

        elif n ==3:
            x = self.fDict['R_eq_s3'](self.K_RC,self.K_CP,self.m_P) /K
            y = 0.
            z = self.fDict['P_eq_s3'](self.K_RC,self.K_CP,self.m_P) /K
        return np.array([x,y,z])
             
    def getInv(self):

        I1 = self.fDict['I_C_s2'](self.K_RC,self.K_CP,self.m_P)
        I2 = self.fDict['I_P_s3'](self.K_RC,self.K_CP,self.m_P)
        I3 = self.fDict['I_P_s4'](self.K_RC,self.K_CP,self.m_P)
        I4 = self.fDict['I_C_s5'](self.K_RC,self.K_CP,self.m_P)

        return I1,I2,I3,I4


    def updateFocalParam(self,focalParam,val):
        if focalParam == 'K_CP':
            self.K_CP = val
        elif focalParam == 'K_RC':
            self.K_RC = val
        else:
            self.m_P = val

    def AssemblySimulation(self,type,N = 5):

        # First Invasion
        if type == 'Cfirst':
            P = self.getEqScenario(2)
        else:
            P = self.getEqScenario(3)
        
        # Second Invasion
        self.UpdateInitCondition(P,type)
        self.setinitCondition(P)
        
        run2 = self.Simulate()

        if run2[1]['message'] == 'Excess work done on this call (perhaps wrong Dfun type).':
            
            run2,r = self.Simulate2(N=N)
            if not(r.successful()):
                run2,r = self.Simulate2(N=N,integrator='dopri5')

            run2 = format(run2,N)
        else:
            run2 = format(run2[0],N)

        return run2
        
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


def format(R,N,tol =1e-12):
    for i in xrange(-N,1,1):
        Positiveformat(R[i],tol)
    return np.array(R[-N:])

def Positiveformat(Array,tol):
    for i in xrange(len(Array)):
        if Array[i]< tol :
            Array[i] = 0.
    return Array    

def SaveData(F,direction,header,footer):
    newF=np.array(FormatData(F)).transpose()
    OutputF = OutputInvData(newF,header,footer,'','')
    OutputF.WriteInvasibility(direction,delimiter=',')
    
  
    
def FormatData(F):   
    x = []
    y = []
    z = []
    
    for i in range(len(F)):
        xhold =[]
        yhold =[]
        zhold =[]
        
        for k in range(len(F[i])):
            xhold.append(MyFloat(F[i][k][0]))
            yhold.append(MyFloat(F[i][k][1]))
            zhold.append(MyFloat(F[i][k][2]))
        
        if len(xhold)>1:
            xhold =MyTuple(xhold)
            yhold =MyTuple(yhold)
            zhold =MyTuple(zhold)
    
        x.append(xhold)
        y.append(yhold)
        z.append(zhold)
            

    return [x,y,z]
