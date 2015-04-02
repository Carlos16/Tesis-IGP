from InvasionAnalysis import *



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


    

    
