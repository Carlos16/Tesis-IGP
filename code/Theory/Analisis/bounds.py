#Get the zero points of a function 
from roots import *
#k_exp=symbols('k_exp')


def Get_bounds(f,get_roots,m_range,upper_guess,lower_guess,k_2_range=[],k_sim = True,debug=False,sep=0.05):
    if not k_sim:
        K_CP_d={}
        for mass in m_range:
            if debug:
                print mass
            K_CP_d[mass]=[]
            for k_2 in k_2_range:
                if debug:
                    print k_2
                new_f = lambdify(k_exp,f(10**k_exp,k_2,mass))
                K_CP_d[mass].append(get_roots(new_f,range_guess,debug=debug))
        return K_CP_d
    else:
        K_CP_l=[]
        for i in range(len( m_range)):
            K_CP_l.append(Roots(m_range[i],lower_guess[i],upper_guess[i],sep,f,get_roots))
        return K_CP_l




def Roots(mass,lower_guess,upper_guess,sep,f,get_roots):
    range_guess = np.arange(lower_guess,upper_guess,sep)

    return get_roots(f,mass,10**(range_guess))





def Get_bounds2(f,get_roots,x_range,search_range,additionalPar =0, k_sim = True,debug=False):
    """
    Search for the zeros of the function f  using the brentq optimization algorithm, separates two cases: \n
    * k_sim is true : f is a bivariate function and the additional argument is given by x_range[i]
    * k_sim is false: f is a trivariate function and the additinal arguments are given by (x_range[i],additionalPar)
    """
    Zeros = []
    if not k_sim:
        for i in range(len(x_range)):
            args = [x_range[i],additionalPar]
            Zeros.append(get_roots(f,args,search_range))              
    else:
        for i in range(len(x_range)):
            args = [x_range[i]]
            Zeros.append(get_roots(f,args,search_range))
    
    return Zeros

### code for plotting and handling the output assuming that each mass point is visited only once if
### we view the boundary as a function of time 
def procce_(Arr,Arr2):
    Y=[]
    X=[]
    for i in range(len(Arr)):
        n = len(Arr[i])
        for k in range(n):
            try:
                Y[k].append(Arr[i][k])
                X[k].append(Arr2[i])
            except:
                X.append([])
                Y.append([])
                Y[k].append(Arr[i][k])
                X[k].append(Arr2[i])
            
    return X,Y

def plot_(x,y,ax,label,color='b',log=True,marker='o',linestyle='-'):
    n=len(x)
    for i in range(n):
        if log:
            if i==n-1:
                ax.plot(np.log10(x[i]),y[i],color=color,marker=marker,linestyle=linestyle,label=label)
            ax.plot(np.log10(x[i]),y[i],color=color,marker=marker,linestyle=linestyle)
    

class Cubic(object):
    def __init__(self,a,b,c,d):
        self.d = d
        self.c = c
        self.b = b
        self.a = a
        self.Xn = -float(b)/(3*a)
        self.Yn = 1
        self.delta2 = 1
        self.l = 1
        self.h2 = 1
        self.delta = 1
        self.h = 1
        self.geometricdiscriminant = 0.
        
    def Evaluate(self,x):
        return self.a*x**3 + self.b*x**2 + self.c*x+self.d
    
    def derivativeEvaluate(self,x):
        return 3*self.a*x**2 + 2*self.b*x + self.c 
    def setYn(self):
        #self.Yn = 2*self.b**3/(27*self.a**2) - self.b*self.c/(3*self.a) + self.d
        self.Yn = self.Evaluate(self.Xn)
    
    def getGeometricDiscriminant(self):
        self.geometricdiscriminant = self.Yn**2 - self.h2
    def setdelta2(self):
        self.delta2 = (self.b**2 - 3.*self.a*self.c)/(9.*self.a**2)
    
    def seth2(self):
        self.h2 = 4.*self.a**2*self.delta2**3
        
    def getRealRoots(self):
        GD = self.geometricdiscriminant
        if GD >0 :
            A1 = 1./(2*self.a)*(-self.Yn + np.sqrt(GD))
            A2 = 1./(2*self.a)*(-self.Yn - np.sqrt(GD))
            return [self.Xn + np.sign(A1)*np.abs(A1)**(1./3) +np.sign(A2)* np.abs(A2)**(1./3)]
        
        if GD == 0:
            if self.Yn == 0:
                return self.Xn
            else:
                A = self.Yn/(2*self.a)
                self.delta = np.sign(A)*np.abs(A)**(1./3)
                
                return [self.Xn + self.delta, self.Xn - 2*self.delta]
        else:
            self.delta =np.sign(self.a)*np.sqrt(self.delta2)
                
            self.h = np.sqrt(self.h2)
            
            theta = np.arccos(-self.Yn/self.h)/3.
            
            alfa = self.Xn + 2*self.delta*np.cos(theta)
            beta = self.Xn + 2*self.delta*np.cos(theta+2*np.pi/3)
            gamma = self.Xn + 2*self.delta*np.cos(theta+4*np.pi/3)
            
            return [alfa,beta,gamma]
        
    


