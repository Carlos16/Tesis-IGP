from SimulationDynamics import *


params1={'w':0.75,'pd':0.21,'pv':0.26,'Er':0.,'Ek':0.,'ER':0.,
        'EC':0.,'EP':0.,'Eq1':0.,'Eq2':0.,'TR':300.,'TC':300,
        'TP':300.,'D_R':2.,'D_C':2.,'fmC':'Grazing','thermyR':'Ecto',
        'thermyC':'Ecto','thermyP':'Ecto','fmPC':'Active','fmPR':'Grazing','k0':1.,'r0':1.71e-6,'a012':10**(-3.08),
        'a03':10**(-3.08),'d0':3.,'q10':4.15e-8,'q20':4.15e-8,'v0R':1.,'v0C':1,'v0P':1,'k':b_k,'e1':0.3,'e2':0.1,'e3':0.5,
        'hC0':1.,'hP0':1.,'formPC':2,'formPR':2,'formC':2,'a':1.,'b':0.02,'c':0}

params2={'w':0.75,'pd':0.2,'pv':0.26,'Er':0.,'Ek':0.,'ER':0.,
        'EC':0.,'EP':0.,'Eq1':0.,'Eq2':0.,'TR':300.,'TC':300,
        'TP':300.,'D_R':3.,'D_C':3.,'fmC':'Grazing','thermyR':'Ecto',
         'thermyC':'Ecto','thermyP':'Ecto','fmPC':'Active','fmPR':'Grazing','k0':300.,'r0':1.71e-6,'a012':10**(-1.77),
         'a03':10**(-1.77),'d0':3.,'q10':4.15e-8,'q20':4.15e-8,'v0R':1.,'v0C':1,'v0P':1,'k':b_k,'e1':0.3,'e2':0.1,'e3':0.5,
        'hC0':1.,'hP0':1.,'formPC':2,'formPR':2,'formC':2,'a':1.,'b':0.02,'c':0}


def Test(params,ksim,xRange,K_RC,K_CP,m_P,focalParam,ParamRange,AssemblyType,type='LV'):
    Prueba = BSR(params,type,xRange,ksim)
    Prueba.setfDict()
    Assembly = Dynamics(Prueba,0,1e5,1,K_RC,K_CP,m_P)
    
    Bif = Assembly.Bifurcation(focalParam,ParamRange,AssemblyType)
    
    return Bif



if __name__== '__main__':
    Test(params1,False,[-10,5],2e0,1e3,1e-5,'K_CP',10**np.arange(-10,5,0.1),type='LV')

    
    
    
    
    
    
    

