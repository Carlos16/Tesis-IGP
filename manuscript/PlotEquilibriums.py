from Plotf1 import *

w = 0.75
pd1 = 0.2
pd2 =0.21
pv = 0.26
e1 = 0.3
e2 = 0.1
e3 = 0.5
k01 = 30
k02 = 0.1
alfa01 = 10**(-1.77)
alfa02 = 10**(-3.08)
b1 = 1.5
KRC =10**(np.arange(-13,6,0.01))
mC = 1e-3
q10 = 4.15e-8
q20 = 4.15e-8
r0 = 1.71e-6

fc = 'Grazing'


h1 = pv + 4*pd1 
h2 = pv + 2*pd2 



f11 = g(KRC,pv,pd1,1.,1.,0.,0.,3,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b1,0)
f12 = g(KRC,pv,pd2,1.,1.,0.,0.,2,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b1,0)
d01 = r0/(alfa01 * KRC**(1-w) *f11)
d02 = r0/(alfa02 * KRC**(1-w) *f12)
d11 = q10/(e1*alfa01*k01*KRC**(1-w)*f11)
d12 = q10/(e1*alfa02*k02*KRC**(1-w)*f12)

g11 = g(KRC,pv,pd1,1.,1.,0.,0.,3,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b1,0)* KRC**(1-w)
g12 = g(KRC,pv,pd2,1.,1.,0.,0.,2,1.,1.,fc,'Ecto','Ecto',1.)*f(KRC,2,1,b1,0)* KRC**(1-w)

Ceq1 = d01*mC**(w-h1)*(1 - d11*mC**(2*w - h1 -1))
Ceq2 = d02*mC**(w-h2)*(1 - d12*mC**(2*w - h2 -1))

Req1 = q10*mC**(w - h1)/(e1*alfa01*f11) 
Req2 = q10*mC**(w - h2)/(e1*alfa02*f12) 


#Max = (((1 + 2*h1 - 3*w)*d11)/(h1-w))**(1/(1 + h1 - 2*w))

fig,axes = plt.subplots(ncols=2)

mpl.rcParams['text.usetex']= True



axes[0].plot(np.log10(KRC) , g11)
axes[0].plot(np.log10(KRC) , g12)

"""
axes[0].plot(np.log10(KRC) , np.log10(Req1))
axes[0].plot(np.log10(KRC) , np.log10(Req2))
axes[0].set_xlim([-10,6])
axes[0].set_xticks([-10,6])
axes[0].set_xlabel(r'$\log_{10}(k_{RC})$')
axes[0].set_ylabel(r'$\log_{10}(R_{Eq})$')
"""
axes[1].plot(np.log10(KRC) , np.log10(Ceq1))
axes[1].plot(np.log10(KRC) , np.log10(Ceq2))
#axes[1].plot([np.log10(Max),np.log10(Max)],[-4,1],'--')
#axes[1].set_ylim([-4,1])
axes[1].set_xlabel(r'$\log_{10}(k_{RC})$')
axes[1].set_ylabel(r'$\log_{10}(C_{Eq})$')
axes[1].set_xticks([-10,6])
axes[1].set_yticks([-3,0,1])
#axes[1].set_xticklabels(["-10",r'$m_{max}$','6'])

plt.show()
#plt.savefig('./Plots/RCeqKRC.pdf')



