import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def F(p,h):
    return ( h/(p - h))**(1/p)


pd1 = 0.2
pd2 = 0.21
p1 = np.concatenate([np.arange(0.42,2,0.01),np.arange(2,100,0.1)])
p2 = np.concatenate([np.arange(0.235,2,0.01),np.arange(2,100,0.1)])

fig,axes = plt.subplots(ncols=1)


mpl.rcParams['text.usetex']= True
axes.plot(np.log10(p1),np.log10(F(p1,2*pd1)),'b-')
axes.plot(np.log10(p2),np.log10(F(p2,pd2)),'r-')
axes.plot([np.log10(2*pd1),np.log10(2*pd1)],[-1,3],'b--')
axes.plot([np.log10(pd2),np.log10(pd2)],[-1,3],'r--')
axes.set_ylim([-1,3])
axes.set_xlabel(r'$\log_{10}(\phi)$')
axes.set_ylabel(r'$\log_{10}(k^*)$')
axes.set_xticks([-1.0,np.log10(pd2),np.log10(2*pd1),0.5,2.0])
axes.set_xticklabels(["-1.0",r"$h_{D_2}$",r"$h_{D_3}$","0.5","2.0"])
#plt.show()
plt.savefig("kmaxGrazing.pdf")
