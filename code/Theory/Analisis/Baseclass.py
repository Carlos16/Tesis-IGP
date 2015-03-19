from functions import *
from Auxiliaryfunctions import *
from OutputClasses import *
import RMEquibrium as RME


class BSR(object):
    """
    Class that stores the basic information for the analysis such as:
    * The mode of the ODE system(RM or LV)
    * The approach(Top down or bottom up)
    * Generate the dictionary with all the functions and the x and y ranges that will be explored.
    * ksim is a boolean variable that set the assumption about similar or not body size ratios across trophic levels.
    * In the case of the Active-Grazing-Grazing strategy, it exists a singular point in the x axis when ksim is false, in this case we set a finer xrange
    in points near to it, this is stored in the xFocus_sep variable.
    """
    def __init__(self,params,mode,xLims, ksim = True):
        self.params = params
        self.mode = mode
        self.bottom = False
        self.xLims = xLims
        self.xFocus = [0,0]
        self.NormalX_sep = 0.01
        self.xFocus_sep = 0.01
        self.Y_sep = 0.05
        self.YLims = [-10,5]
        self.fDict = {}
        self.ksim = ksim
        
        self.K_CP = symbols('K_CP')
        self.K_RC = symbols('K_RC')
        self.m_P = symbols('m_P')
        self.C = symbols('C')
        self.R = symbols('R')
        self.P = symbols('P')
            
    def setfDict(self):
        """
        Constructs the dictionary containing all the elementary functions used in the analysis
        Separates two different approachs:
        if bottom is true  m_C and m_P are expressed in terms of m_R, else m_C and m_R are expressed in terms of m_P, in reality the scenarios will give
        similar pictures, we are just rotating the plane if we interchange them        
        """
        bottom = self.bottom
        if bottom:
            self.fDict = func_transform(self.params,self.K_CP,self.K_RC,self.m_R/(self.K_RC*self.K_CP),self.R,self.C,self.P
                                        ,bottom = bottom,sim=self.ksim)
        else:
            self.fDict  = func_transform(self.params,self.K_CP,self.K_RC,self.m_P,self.R,self.C,self.P,sim=self.ksim)
    def setInvFunctions(self):
        mode = self.mode
        if mode == "RM":
            return ['I_C_s2RM','I_P_s3RM','I_P_s4RM','I_C_s5RM']
        else:
            return ['I_C_s2','I_P_s3','I_P_s4','I_C_s5']
    def setxRange(self):
        """
        Sets the range of X points which are going to be explored
        """
        A = np.arange(self.xLims[0],self.xFocus[0],self.NormalX_sep)
        B = np.arange(self.xFocus[0],self.xFocus[1],self.xFocus_sep)
        C = np.arange(self.xFocus[1],self.xLims[1],self.NormalX_sep)
        D = np.concatenate([A,B,C])
        return 10**(D)
    
    def getandSetxFocus(self,nPoints,dist):
        """
        nPoints : the number of points in which the interval [x0-dist,x0] is going to be divided, x0 is where the singularity occurres (if it exists)
        """
        e1 = self.params['e_1']
        e2 = self.params['e_2']
        q10 = self.params['q10']
        q20 = self.params['q20']
        pd = self.params['pd']
        pv = self.params['pv']
        D_R = self.params['D_R']
        w = self.params['w']
        xSing = np.log10((e2*q10/(e1*q20))**(1/(pd*(D_R-1) + pv - w)))
        lowEnd = xSing-dist
        self.setxFocus([lowEnd,xSing])
        self.setxFocusSep(dist/nPoints)

    def setxFocus(self,focusLims):
        self.xFocus = focusLims
    def setxFocusSep(self,Focus_sep):
        self.xFocus_sep = Focus_sep
    def getParams(self):
        return self.params
            
    def getfDict(self):
        """
        return the dictionary of functions
        """
        return self.fDict

    def getmode(self):
        """
        return the mode of the analysis, possible answers at the moment are LV and RM
        """
        return self.mode
    def getxLims(self):
        """
        return the boundary points of the xRange list(assumed as an Interval)
        """
        return self.xLims
    
    def setyRange(self):
        Y_logrange= np.arange(self.YLims[0],self.YLims[1],self.Y_sep)
        return 10**Y_logrange
        
