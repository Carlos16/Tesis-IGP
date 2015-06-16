import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
def g(k_,pv,pd,T1,T2,E1,E2,D,v01,v02,fm,thermy1,thermy2,k):
    if thermy1=="Ecto" and thermy2=="Ecto":
        if fm=="Active":
            delta=np.exp(-(1./k)*(E1/T1 - E2/T2))
            return v02*np.exp(-E2/(k*T2))*(1+(v01/v02)*k_**(2*pv)*delta**2)**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)*np.exp(-E1/(k*T1))
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)*np.exp(-E2/(k*T2))
    elif thermy1=="Endo" and thermy2=="Ecto":
        if fm=="Active":
            return v02*(1+(v01/v02)*k_**(2*pv)*np.exp(-2*E1/(k*T1)))**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)*np.exp(-E1/(k*T1))
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)
    elif thermy1=="Ecto" and thermy2=="Endo":
        if fm=="Active":
            return v02*(1+(v01/v02)*k_**(2*pv)*np.exp(2*E2/(k*T2)))**0.5*k_**((D-1)*pd)
        elif fm=="Sit":
            return v01*k_**(pv+(D-1)*pd)
        elif fm=="Grazing":
            return v02*k_**((D-1)*pd)*np.exp(-E2/(k*T2))

def f(k_,form,a,b,c):
    if form == 1:
        return a
    elif form == 2:
        return a/(1+k_**b)
    elif form == 3 :
        return a*exp(-(log(k_)-b)**2/(2*c))




pd1 = 0.21
pd2 = 0.2
pv = 0.26
b1 = 0.1
b2 = 0.5
b3 = 0.3
KRC =10**(np.arange(-13,6,0.01))

fc = 'Grazing'

fig,axes = plt.subplots(ncols=2,sharex = True)
#mpl.rcParams['text.usetex']= True

axes[0].plot(np.log10(KRC),g(KRC,pv,pd1,1.,1.,0.,0.,3,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b1,0)*KRC**(0.25),'b-')
axes[0].plot(np.log10(KRC),g(KRC,pv,pd2,1.,1.,0.,0.,2,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b1,0)*KRC**(0.25),'r-')

axes[0].set_ylim([0,10])

axes[0].set_xlabel(r'$\log_{10}(k_{ij})$')
axes[0].set_ylabel(r'$f$')

axes[1].plot(np.log10(KRC),g(KRC,pv,pd1,1.,1.,0.,0.,3,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b2,0)*KRC**(0.25),'b')
axes[1].plot(np.log10(KRC),g(KRC,pv,pd2,1.,1.,0.,0.,2,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b3,0)*KRC**(0.25),'r')

axes[1].set_xlabel(r'$\log_{10}(k_{ij})$')
axes[0].set_xlim([-10,6])

plt.show()
#plt.savefig('./Plots/f1Sit.pdf')
