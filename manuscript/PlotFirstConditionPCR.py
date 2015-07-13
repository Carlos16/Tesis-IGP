from InvasibilityBFunctions import *

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



#3D

I51 = InvP4B(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,0.02,D_R,w,'Grazing','Active','Grazing')
I61 = InvP4B(KRC,KCP,e11,e21,e3,alfa0,pv,pd,k0,r0,q10,q20,0.02,D_R,w,'Sit','Sit','Active')
I71 = InvP4B(KRC,KCP,e11,e21,e3,alfa0,pv,pd,k0,r0,q10,q20,0.02,D_R,w,'Active','Active','Active')


I52 = InvP4B(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,0.2,D_R,w,'Grazing','Active','Grazing')
I62 = InvP4B(KRC,KCP,e11,e21,e3,alfa0,pv,pd,k0,r0,q10,q20,0.2,D_R,w,'Sit','Sit','Active')
I72 = InvP4B(KRC,KCP,e11,e21,e3,alfa0,pv,pd,k0,r0,q10,q20,0.2,D_R,w,'Active','Active','Active')

I53 = InvP4B(KRC,KCP,e1,e2,e3,alfa0,pv,pd,k0,r0,q10,q20,2,D_R,w,'Grazing','Active','Grazing')
I63 = InvP4B(KRC,KCP,e11,e21,e3,alfa0,pv,pd,k0,r0,q10,q20,2,D_R,w,'Sit','Sit','Active')
I73 = InvP4B(KRC,KCP,e11,e21,e3,alfa0,pv,pd,k0,r0,q10,q20,2,D_R,w,'Active','Active','Active')


#2D

I512 = InvP4B(KRC,KCP,e1,e2,e3,alfa02,pv,pd2,k02,r0,q10,q20,0.02,D_R2,w,'Grazing','Active','Grazing')
I612 = InvP4B(KRC,KCP,e11,e21,e3,alfa02,pv,pd2,k02,r0,q10,q20,0.02,D_R2,w,'Sit','Sit','Active')
I712 = InvP4B(KRC,KCP,e11,e21,e3,alfa02,pv,pd2,k02,r0,q10,q20,0.02,D_R2,w,'Active','Active','Active')


I522 = InvP4B(KRC,KCP,e1,e2,e3,alfa02,pv,pd2,k02,r0,q10,q20,0.2,D_R2,w,'Grazing','Active','Grazing')
I622 = InvP4B(KRC,KCP,e11,e21,e3,alfa02,pv,pd2,k02,r0,q10,q20,0.2,D_R2,w,'Sit','Sit','Active')
I722 = InvP4B(KRC,KCP,e11,e21,e3,alfa02,pv,pd2,k02,r0,q10,q20,0.2,D_R2,w,'Active','Active','Active')

I532 = InvP4B(KRC,KCP,e1,e2,e3,alfa02,pv,pd2,k02,r0,q10,q20,2,D_R2,w,'Grazing','Active','Grazing')
I632 = InvP4B(KRC,KCP,e11,e21,e3,alfa02,pv,pd2,k02,r0,q10,q20,2,D_R2,w,'Sit','Sit','Active')
I732 = InvP4B(KRC,KCP,e11,e21,e3,alfa02,pv,pd2,k02,r0,q10,q20,2,D_R2,w,'Active','Active','Active')



from matplotlib import cm

fig,axes = plt.subplots(ncols=3,nrows=2,sharex = True,sharey =True)
levs= [-1,0]




axes[0,0].contour(np.log10(KCP),np.log10(KRC),I512,levels=levs,cmap = cm.Greens)
axes[0,0].contour(np.log10(KCP),np.log10(KRC),I612,levels=levs,cmap = cm.Reds)
axes[0,0].contour(np.log10(KCP),np.log10(KRC),I712,levels=levs,cmap =cm.Blues)

axes[0,1].contour(np.log10(KCP),np.log10(KRC),I522,levels=levs,cmap = cm.Greens)
axes[0,1].contour(np.log10(KCP),np.log10(KRC),I622,levels=levs,cmap = cm.Reds)
axes[0,1].contour(np.log10(KCP),np.log10(KRC),I722,levels=levs,cmap =cm.Blues)

axes[0,2].contour(np.log10(KCP),np.log10(KRC),I532,levels=levs,cmap = cm.Greens)
axes[0,2].contour(np.log10(KCP),np.log10(KRC),I632,levels=levs,cmap = cm.Reds)
axes[0,2].contour(np.log10(KCP),np.log10(KRC),I732,levels=levs,cmap =cm.Blues)


axes[1,0].contour(np.log10(KCP),np.log10(KRC),I51,levels=levs,cmap = cm.Greens)
axes[1,0].contour(np.log10(KCP),np.log10(KRC),I61,levels=levs,cmap = cm.Reds)
axes[1,0].contour(np.log10(KCP),np.log10(KRC),I71,levels=levs,cmap =cm.Blues)

axes[1,1].contour(np.log10(KCP),np.log10(KRC),I52,levels=levs,cmap = cm.Greens)
axes[1,1].contour(np.log10(KCP),np.log10(KRC),I62,levels=levs,cmap = cm.Reds)
axes[1,1].contour(np.log10(KCP),np.log10(KRC),I72,levels=levs,cmap =cm.Blues)

axes[1,2].contour(np.log10(KCP),np.log10(KRC),I53,levels=levs,cmap = cm.Greens)
axes[1,2].contour(np.log10(KCP),np.log10(KRC),I63,levels=levs,cmap = cm.Reds)
axes[1,2].contour(np.log10(KCP),np.log10(KRC),I73,levels=levs,cmap =cm.Blues)

#plt.show()
axes[1,1].set_xlabel(r'$\log_{10}(k_{CP})$')
axes[0,0].set_ylabel(r'$\log_{10}(k_{RC})$')
axes[1,0].set_ylabel(r'$\log_{10}(k_{RC})$')
plt.savefig('NecessityPCR.pdf')





