from InvasionAnalysis import *

params = {'formC': 2, 'fbs': 'comb1', 'TR': 300.0, 'fmC': 'Grazing', 'd0': 3.0, 'e3': 0.5, 'e2': 0.1, 'e1': 0.3, 'pv': 0.26, 'thermyC': 'Ecto', 'D_R': 2.0, 'fm': 'comb1', 'fmPR': 'Grazing', 'a03': 0.0008317637711026709, 'thermyR': 'Ecto', 'pd': 0.21, 'D_C': 2.0, 'thermyP': 'Ecto', 'Eq1': 0.0, 'Eq2': 0.0, 'e': 'comb1', 'b':0.02, 'formPR': 2, 'formPC': 2, 'q20': 4.15e-08, 'a012': 0.0008317637711026709, 'D': '2D', 'v0P': 1, 'v0R': 1.0, 'hP0': 1.0, 'TP': 300.0, 'k0': 0.01, 'v0C': 1, 'TC': 300, 'hC0': 1.0, 'fmPC': 'Active', 'EC': 0.0, 'q10': 4.15e-08, 'EP': 0.0, 'ER': 0.0, 'a': 1.0, 'c': 0, 'r0': 1.71e-06, 'Ek': 0.0, 'k': 8.617332400000001e-05, 'w': 0.75, 'Er': 0.0}

DirInv = 'c:/Users/Carlos/Documents/Thesis/Tesis/modules/InvPrueba.csv'
DirZB = 'c:/Users/Carlos/Documents/Thesis/Tesis/modules/ZBPrueba.csv'
DirW = 'c:/Users/Carlos/Documents/Thesis/Tesis/modules/WidthsPrueba.csv'

def Test(params,ksim,xRange,DirInv,DirZB,DirW,type='LV',mass=1e-10):
    Prueba = BSR(params,type,xRange,ksim)
    Prueba.setfDict()
    HeaderInv = ['Inv C2','Inv P3','Inv P4','Inv C5']
    Inv  = InvBoundaries(Prueba,mass)
    Inv.setAndWriteInvBoundaries(HeaderInv,DirInv) 
    Inv.setPositiveBoundaries()
    Inv.setAndWriteWidthsZones(DirW,DirZB)


if __name__== '__main__':
    Test(params,False,[-10,5],DirInv,DirZB,DirW,type='LV')


    

    

