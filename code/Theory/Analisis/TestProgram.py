from InvasionAnalysis import *
params1={'w':0.75,'pd':0.21,'pv':0.26,'Er':0.,'Ek':0.,'ER':0.,
        'EC':0.,'EP':0.,'Eq1':0.,'Eq2':0.,'TR':300.,'TC':300,
        'TP':300.,'D_R':2.,'D_C':2.,'fmC':'Grazing','thermyR':'Ecto',
        'thermyC':'Ecto','thermyP':'Ecto','fmPC':'Active','fmPR':'Grazing','k0':0.01,'r0':1.71e-6,'a012':10**(-3.08),
        'a03':10**(-3.08),'d0':3.,'q10':4.15e-8,'q20':4.15e-8,'v0R':1.,'v0C':1,'v0P':1,'k':b_k,'e1':0.3,'e2':0.1,'e3':0.5,
        'hC0':1.,'hP0':1.,'formPC':2,'formPR':1,'formC':1,'a':1.,'b':2.,'c':0}

params2={'w':0.75,'pd':0.2,'pv':0.26,'Er':0.,'Ek':0.,'ER':0.,
        'EC':0.,'EP':0.,'Eq1':0.,'Eq2':0.,'TR':300.,'TC':300,
        'TP':300.,'D_R':3.,'D_C':3.,'fmC':'Grazing','thermyR':'Ecto',
         'thermyC':'Ecto','thermyP':'Ecto','fmPC':'Active','fmPR':'Grazing','k0':3.,'r0':1.71e-6,'a012':10**(-1.77),
         'a03':10**(-1.77),'d0':3.,'q10':4.15e-8,'q20':4.15e-8,'v0R':1.,'v0C':1,'v0P':1,'k':b_k,'e1':0.3,'e2':0.1,'e3':0.5,
        'hC0':1.,'hP0':1.,'formPC':2,'formPR':1,'formC':1,'a':1.,'b':2.,'c':0}

DirInv = 'c:/Users/Carlos/Documents/Thesis/Tesis/modules/InvPrueba.csv'
DirZB = 'c:/Users/Carlos/Documents/Thesis/Tesis/modules/ZBPrueba.csv'
DirW = 'c:/Users/Carlos/Documents/Thesis/Tesis/modules/WidthsPrueba.csv'

def Test(params,ksim,xRange,DirInv,DirZB,DirW,type='LV',mass=1e+05):
    Prueba = BSR(params,type,xRange,ksim)
    Prueba.setfDict()
    HeaderInv = ['Inv C2','Inv P3','Inv P4','Inv C5']
    Inv  = InvBoundaries(Prueba,mass)
    Inv.setAndWriteInvBoundaries(HeaderInv,DirInv) 
    Inv.setPositiveBoundaries()
    Inv.setAndWriteWidthsZones(DirW,DirZB)


if __name__== '__main__':
    Test(params,False,[-10,5],DirInv,DirZB,DirW,type='LV')


    

    
