from InvasibilityBFunctions import *
from numpy import ma
KCP = 10**np.arange(-13,10,0.01)
KRC =10**(np.arange(-13,10,0.01))
KCP,KRC = np.meshgrid(KCP,KRC)

#Dimension independent parameters
w = 0.75
b = 0.02
q10 = 4.15e-8
q20 = 4.15e-8
r0 = 1.71e-6
pv = 0.26


#3D
k0 = 30.
alfa0 = 10**(-1.77)
D_R = 3.
pd = 0.2


#2D
D_R2 = 2.
k02 = 0.1
alfa02 = 10**(-3.08)
pd2 = 0.21




#Gr-Gr-Ac
e1 = 0.3
e2 = 0.1
e3 = 0.5


#Other

e11 = 0.6
e21 = 0.4


def PlotBoundsDistancesI5(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):
    
    I1 = InvP5(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc)
    I2 = InvP3(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc)

    return np.exp(I1)


def PlotBoundsDistancesCoex(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc):

    
    
    I1 = InvP5(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc)
    I2 = InvP4(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,b,D,w,fpr,fpc,fc)

    return np.exp(I1) - np.exp(I2)


fig,axes = plt.subplots()

I2 = PlotBoundsDistancesCoex(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,2.,D_R,w,'Grazing','Active','Grazing')
"""
I1 = InvP5(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,2.,D_R,w,'Grazing','Active','Grazing')

I3 = InvP4(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,2.,D_R,w,'Grazing','Active','Grazing')
#I =  ma.masked_array(I,np.isnan(I))


#print(InvP5(10**-7.5,10**9.5,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,2.,D_R,w,'Grazing','Active','Grazing'))



c = axes.contourf(np.log10(KCP),np.log10(KRC),np.exp(I1),norm = m.colors.LogNorm())
cbar = plt.colorbar(c)
fig2,axes2 = plt.subplots()



fig3,axes3 = plt.subplots()
c3 = axes3.contourf(np.log10(KCP),np.log10(KRC),np.exp(I3),norm = m.colors.LogNorm())
cbar = plt.colorbar(c3)
"""
c = axes.contourf(np.log10(KCP),np.log10(KRC),I2,norm = m.colors.LogNorm())
cbar = plt.colorbar(c)


plt.show()


