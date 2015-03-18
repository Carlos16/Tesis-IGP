from functions import *
from MTP import *
from neweq import *
from scipy.integrate import odeint
from scipy.optimize import fsolve
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
        
class InvBoundaries(BSR):
    """ 
    Class that stores and computes the values for the invasibility boundaries of the distinct scenarios, it receives as initial
    input the values of the parameters and the invasibility functions 
    """
    def __init__(self,workingData,currentMass= 0.):
        
        BSR.__init__(self,workingData.getParams(),workingData.getmode(),workingData.getxLims(),workingData.ksim)
        self.xRange = workingData.setxRange()
        self.yRange = workingData.setyRange()
        self.UpGuess =  10
        self.LowGuess = -10
        self.guessSep = 0.03
        self.InvFunctions=workingData.setInvFunctions()
        self.InvBounds={}
        self.UnEditedInvBounds={}
        self.fDict = workingData.getfDict()
        self.EqValues = []
        self.stabVals = []
        self.MTPvals = []
        self.distance = 1e-12
        self.PositiveBoundaries={}
        self.currentMass = currentMass
        self.xFocus = workingData.xFocus
        self.xFocus_sep = workingData.xFocus_sep
        self.DBound = 0
        self.EfDif = self.params['e_1']*self.params['e_3'] - self.params['e_2']

        self.Widths={}
        self.SZ={}
        self.SZBounds={}
        self.Footer = constructFooter(self.params)
            
    
    def setUpGuess(self,Guess):
        self.UpGuess = Guess
            
    def UpdateMass(newMass):
        self.currentMass = newMass
        
    def setAndWriteInvBoundaries(self,Header,Direction):
        """
        * Find the zero boundaries for each of the functions specified in self.InvFunctions
        * Write them to a csv file specified in Direction.
        
        Header refers to the first row of the csv file.
        """

        self.setInvBoundaries()
        self.writeInvasibilityValues(Header,Direction)
        
        
    def setInvBoundaries(self):
        """ 
        For each of the functions present in the InvFunctions list, computes the invasibility boundaries by means of the
        Scipy.Brentq numerical method.
        For each of the invasibility functions it returns a dict object with x and y keys indicating the zeros(size ratio values delimiting
        the boundaries) for the function at each x (mass), the values of x and y are a list of lists. the total number of sublists denotes the maximum
        number of zeros found at any location x.
        if e1e3 - e2 = EfDif <0 it also computes the boundaries for the D function 
        """
        searchRange = 10**(np.arange(self.LowGuess,self.UpGuess,self.guessSep))
        xRange = self.xRange
        if self.EfDif<0:
            self.InvFunctions.append('D')
        for invfunc in self.InvFunctions:
            bound_dict = self.InvBoundary(invfunc,searchRange,xRange)
            self.InvBounds[invfunc] = bound_dict           
    
    def InvBoundary(self,Invfunc,searchRange,xRange):
        """Find the inv boundaries(zeros of the invasibility function) using the brentq method from the SciPy package,
        input arguments, the invasibility function, the limits for the interval to look for zeros(depending on the x values) """
        
        f = self.fDict[Invfunc]
        K = Get_bounds2(f,Get_roots,xRange,searchRange,additionalPar= self.currentMass,k_sim = self.ksim)
        x,y= procce_(K,xRange)
        self.UnEditedInvBounds[Invfunc]=K
        return {'x':x,'y':y}
    
    
    def getBounds(self):
        """ return the dictionary containing all the invasibility boundaries"""
        return self.InvBounds
    
    def writeInvasibilityValues(self,Header,direction,delimiter=','):
        """Write data into a csv file specified in direction
        @param Boundaries a list of 2-elements lists which each of them stores the X and Y coordinates
        of the invasibility boundaries computed in the analysis
        @param the list of params used to compute the boundaries which will be used in the footer of the csv file
        @direction the system direction where the file is going to be stored
        @Header the first row of the csv file """
       
        new_Boundaries,dists = FormatZones(self.InvFunctions,self.InvBounds)
       
        Footer = self.Footer
        if self.EfDif<0:
            T = Header + ['D']
            OutputFile = OutputInvData(new_Boundaries,T,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],dists)
        else:
            OutputFile = OutputInvData(new_Boundaries,Header,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],dists)
        OutputFile.WriteInvasibility(direction,delimiter)

       
    
    def setPositiveBoundaries(self):
        """ Creates a dict storying the boundary points of the set in which each of the criterions is satisfied """
        for Inv in self.UnEditedInvBounds.keys():
            M,P = RME.GetGuessBounds(self.UnEditedInvBounds[Inv],self.xRange,self.fDict[Inv],self.UpGuess,self.LowGuess,self.ksim,self.currentMass)
            
            self.PositiveBoundaries[Inv] = P
            
                
    def getAndWriteSMTPEq(self,DirectionEq,DirectionStab,DirectionMTP):
        """
        In the case of the LV mode(which is used in the UG thesis) , this method serves to compute the Equilibrium Values, eigenvalues and Maximum Trophic
        positions at each of the (i,j) locations specified by the X and Y 2D arrays.

        * Using the information from the Invasibility boundaries it first computes the boundaries in which each of the Invasibility functions are positive
          If EfDif <0 it also takes into account the boundaries for D 
        * It stores each of the results in 2D arrays (e.g. Zeq)
        * Write each of them to csv files specified in the arguments Directions...
        """
        X,Y = np.meshgrid(self.xRange,self.yRange)
        PPos1 = self.PositiveBoundaries['I_P_s3']
        PPos2 = self.PositiveBoundaries['I_P_s4']
        CPos1 = self.PositiveBoundaries['I_C_s2']
        CPos2 = self.PositiveBoundaries['I_C_s5']
        if self.EfDif < 0 :
            DPos = self.PositiveBoundaries['D']
        else:
            DPos = 0
            
        Zeq,Zstab,ZMTP = self.setStabilityMTPEquilibrium2(X,Y,PPos1,CPos1,PPos2,CPos2,DPos)
        
        self.WriteEquilibriumValues(Zeq,DirectionEq)
        self.WriteEquilibriumValues(Zstab,DirectionStab)
        self.WriteEquilibriumValues(ZMTP,DirectionMTP)
        
        
    
    def setStabilityMTPEquilibrium2(self,X,Y,PPos1,CPos1,PPos2,CPos2,DPos):
        """
        Loop through all the elements of X and Y and computes the equilibrium , eigenvalues and MTP at each (X[i,j],Y[i,j]) position.
        Saves them in the respective Z matrix.        
        Distinguish two cases:
        EfDif >=0:
        * use the getRegion function 
        EfDIf <0:
        * use the getRegion2 function
        Ksim True:
        * AddPar is an empty list
        else:
        * AddPar is [self.currentMass]
        """
        m,n = X.shape
        Zeq = np.zeros(X.shape,dtype=object)
        Zstab=np.zeros(X.shape,dtype=object)
        ZMTP =np.zeros(X.shape,dtype=object)  


        Fargs = [PPos1,CPos1,PPos2,CPos2,Zeq,Zstab,ZMTP]

        # account for differences if similar or nonsimilar size ratios assumption is used, in the second case 
        # add the current mass to the list of arguments
        if self.ksim:
            AddPar = []
        else:
            AddPar = [self.currentMass]

        # account for differences in the getRegion function if D <0
        if self.EfDif >= 0:
            RegionFunct = getRegion
            
        else:
            RegionFunct = getRegion2
            Fargs +=[DPos]

        for y in range(n):
            for x in range(m):
                self.EqStabMTPFunction(x,y,X,Y,AddPar,RegionFunct,*Fargs)

        return Zeq,Zstab,ZMTP

    def EqStabMTPFunction(self,x,y,X,Y,AddPar,getR,PP1,CP1,PP2,CP2,Zeq,Zstab,ZMTP,DP=0):
        """
        * Take points X[x,y] and Y[x,y] and computes the equilibrium, eigenvalues and MTP for them
        * Take an additional argument AddPar which is nonempty if ksim is not True
        * PP1,PP2,CP2,CP1 represent boundaries of the positive regions of each of the invasibility funcionts, if DP is specified it also does the same.
        * Add the computed values to the Zeq,Zstab and ZMTP matrices respectively.
        """
        
        xPoint = X[x,y] ; yPoint = Y[x,y]
        args = [yPoint,xPoint] + AddPar
        r = getR(y,Y[x,y],PP1,CP1,PP2,CP2,DP,self.fDict,args)                       
        zeq,zstab,zMTP = self.SetStabEqMTP2(args,r,self.fDict)
        Zeq[x,y] = FormatEq(zeq) ; Zstab[x,y]=FormatStab(zstab) ; ZMTP[x,y] = MyFloat(zMTP)

    def SetStabEqMTP2(self,args,r,fDict):
        #get eq
        Eq = setEq(args,fDict,r,mode=self.mode)
        #get eigenvalues
        argsEstab = Eq+args
        
        Eigenvals = setEigenvals(argsEstab,fDict,r,mode=self.mode)
        
        #get MTPval
        MTPval = setMTP(args,fDict,r,mode=self.mode)
    
        return Eq,Eigenvals,MTPval
        
           
    def setEquilibriumvalue(self,x,y):
        
        if self.mode == "LV":
            F = self.fDict
            R = new_R(y,y,x,F['R_eq'],F['I_P_s4'],F['I_C_s5'],F['I_C_s2'],F['I_P_s3'],F['R_eq_s2'],F['R_eq_s3'],set_R_2,set_R_3,F['K'])
            C = new_C(y,y,x,F['C_eq'],F['I_P_s4'],F['I_C_s5'],F['I_C_s2'],F['C_eq_s2'],set_C_2)
            P = new_P(y,y,x,F['P_eq'],F['I_P_s4'],F['I_C_s5'],F['I_P_s3'],F['P_eq_s3'],set_P_2)
            
            return CustomizeObjects([R,C,P])
        
        else:
            EqRM  = self.setEqRM(x,y)
            return EqRM
             
    def setEqRM(self,x,y):
        Req = setReqRM(self.fDict,y,x)
        Ceq = [setCeqRM(R,self.fDict,y,x) for R in Req]
        Peq = [setPeqRM(Req[i],Ceq[i],self.fDict,y,x) for i in range(len(Req))]   
        return self.FormatedEq(Req,Ceq,Peq)
        
    def FormatedEq(self,Req,Ceq,Peq):
        EqVals = []
        for i in range(len(Req)):
            EqVals.append(MyInnerTuple([MyFloat(Req[i]),MyFloat(Ceq[i]),MyFloat(Peq[i])]))
        return MyTuple(EqVals)
    
    def WriteEquilibriumValues(self,EqValues,Direction,delimiter=','):
        """
        Write the EqValues stored in matrix format into a file, the direction is specified in Direction and the delimiter default is ','
        which means that the default ouput type is a csv file.
        """
        Footer = constructFooter(self.params)
        OutputFile = OutputInvData(EqValues,"",[Footer],[MyTuple(self.xFocus),self.xFocus_sep],"")
        OutputFile.WriteEquilibrium(Direction,delimiter)

    def WriteWidths(self,Direction):
        """
        Write the Width dictionary to a file specified in Direction.
        """
        Header = self.Widths.keys()
        data = FormatWidths(Header,self.Widths)
        Footer = self.Footer
        OutputFile = OutputInvData(data,Header,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],"")
        OutputFile.WriteWidths(Direction,',')
        

    def setAndWriteWidthsZones(self,DirectionW,DirectionZB,DirectionZones):
        """
        * For each x in xRange and any Invasibility function f calculates the sum of the length of all the intervals which are in the Positive region of f
        * Get the Positive Boundaries for all the target zones in the analysis
        * Write both of the above results to a file whose direction is specified in the arguments DirectionW and DirectionZ.
        """
        self.setWidthsAndZones()
        self.WriteWidths(DirectionW)
        self.WriteZonesBounds(DirectionZB)
        self.WriteZones(DirectionZones)
    
    def setWidthsAndZones(self):
        """
        * Calculates the Boundaries of the Zones described in the Study, which are intersection of the Positive Regions of a subset of the Invasibility functions
        * Calculates the widths of each of the Zones
        * Convert them to the format used for the Invasibility functions and save them in the ZBounds dict
        * if EfDif<0 , calculates the boundary for the stable coexistence zone.
        """

        self.setIntersections()
        self.setZonesBoundaries()
        self.setZones()
        self.setWidths()

    def setIntersections(self):
        r"""
        Find the boundary points of the Intersections of the Positive Regions of some invasibility functions, which delineate the zones expressed in the analysis.
        For example :math:`Z(I_{C4}) := Z(I_{C2}) \cap PR_{4}` where :math:`PR_{4}` is the set for which the Invasibility function I_P_s4 is positive.
        We assume that one dimensional subset of the set are always a finite union of \emph{connected} sets.
        """
        
        I1 = GetIntersection(self.PositiveBoundaries['I_C_s5'],self.PositiveBoundaries['I_P_s3'])
        I2 = GetIntersection(self.PositiveBoundaries['I_C_s2'],self.PositiveBoundaries['I_P_s4'])
        I3 = GetIntersection(I1,I2)
        I4 = GetIntersection(self.PositiveBoundaries['I_C_s5'],self.PositiveBoundaries['I_P_s4'])

       
        self.Intersections['Z(IC5)'] = I1
        self.Intersections['Z(IC4)'] = I2
        self.Intersections['MutualInv'] = I3
        self.Intersections['Coexistence'] = I4
        
        if self.EfDif < 0 :
            I5 = GetIntersection(GetComplement(self.PositiveBoundaries['I_C_s5']),GetComplement(self.PositiveBoundaries['I_P_s4']))
            I6 = GetIntersection(I5,GetComplement(self.PositiveBoundaries['D']))
            I7 = GetIntersection(I4,self.PositiveBoundaries['D'])

            self.Intersections['UnstableCoexistence'] = I6
            self.Intersections['Coexistence'] = GetUnion(I6,I7)
           
    def findWidths(self):
        """
        Find the widths of each of the Zones 
        """
        
        F(self.Widths,getWidths,self.PositiveBoundaries,self.Intersections,self.EfDif)
    def findZonesBoundaries(self):
        """
        Find the boundary of the Zones and format them to the same data structure used in the computation of the Invasibility boundaries
        """
        F(self.ZBounds,GetPositiveBoundaries,self.PositiveBoundaries,self.Intersection,self.EfDif,self.xRange)
    def findZones(self):
        """
        From each of the Boundaries, create a 2D array using xRange and yRange and return a 1 0  array , each 1 located at a position in which
        the point is interior to the Positive Region delimited by the Boundary
        """
        F(self.Zones,GetPositiveRegions,self.PositiveBoundaries,self.Intersection,self.EfDif,self.xRange,self.yRange)

    def WriteZonesBounds(self,Direction):
        """
        Write the SZBounds into a file whose pointer is specified in Direction.
        """
        Header = self.ZBounds.keys()
        data,dist = FormatZones(Header,self.SZBounds)
        Footer = self.Footer
        OutputFile = OutputInvData(data,Header,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],dist)
        OutputFile.WriteInvasibility(Direction,',')

#    def WriteZones(self,Direction):

        

        
########## Auxiliary functions######################        
def GetPositiveRegions(Bounds,xRange,yRange):
    X,Y = np.meshgrid(xRange,yRange)
    m,n = X.shape
    Z = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            Z[i,j] = int(inList(Y[i,j],Bounds[j]))
    return Z

            

def F(D,F1,D1,D2,EfD,*args):
    D['Z(IC2)'] = F1(D1['I_C_s2'],*args)
    D['Z(IC3)'] = F1(D1['I_P_s3'],*args)
    D['Z(IC4)'] = F1(D2['Z(IC4)'],*args)
    D['Z(IC5)'] = F1(D2['Z(IC5)'],*args)
    D['MutualInv'] = F1(D2['MutualInv'],*args)
    D['Coexistence'] = F1(D2['Coexistence'],*args)
    
    if EfD < 0:
        D['UnstableCoexistence'] = F1(D2['UnstableCoexistence'],*args)


def GetPositiveBoundaries(PosPoints,xRange):
    """ 
    Input a range of Points and for each x in Points a set of tuples containing 
    the points of the boundary for the set in which the function has positive
    values.
    Output a dictionary in the format of the invasibility boundary dicts 
    """
    Up,Low,X = RME.ExtractPoints(xRange,PosPoints)
    
    newX =[]
    newY =[]
    for i in range(len(X)):
        newX.append(X[i])
        newX.append(X[i])
        newY.append(Low[i])
        newY.append(Up[i])
    return {'x':newX,'y':newY}

              
def FormatWidths(Keys,WidthDict):
    """
    Format Widths using a custom tupple class to get 20 digits printing
    """
    Data = []
    for i in range(len(WidthDict[Keys[0]])):
        row =[]
        for k in Keys:
            row.append(MyTuple([WidthDict[k][i],0]))
        Data.append(row)
    return Data
 
       
def GetComplement(Boundaries,Min,Max):
    r"""
    From a set of Boundary points of a set, return the Boundary Points of its complement. since each one dimensional slide will be a set of intervals it pieces together the complements of each contiguous pair of intervals in each list :math:`((a,b) \cup (c,d))^c = [b,c]`
    """
    ComplementBoundary = [] 
    for boundPoints in Boundaries:
        ComplementBoundary.append(getCompPoints(boundPoints,Min,Max))
    return ComplementBoundary

def getCompPoints(boundPoints,Min,Max):
    if len(boundPoints[0]) == 1:
        return [(Min,Max)]
    else:
        CompPoints = []
        #Little formating
        NewB = boundPoints
        if boundPoints[0][0]!= Min:
            NewB = [(Min,Min)]+NewB
        if boundPoints[-1][1]!= Max:
            NewB+=[(Max,Max)]
            
        for i in range(len(newB)-1):
            CompPoints.append((newB[i][1],newB[i+1][0]))

        if len(CompPoints)==0:
            CompPoints.append((0,))

        return CompPoints
            
            
            
            
        

                
        

def FormatZones(Header,SZBounds):
    """
    Convert each element of the SZBounds into the format used for the InvFunctions dict
    """
    new_Boundaries = []
    dists = []
    for zone in Header:
        bound = SZBounds[zone]
        new_Bound , distB = Convert(bound['x'],bound['y'])
        new_Boundaries.append(new_Bound)
        dists.append(distB)
    new_Boundaries =ConstructArray(new_Boundaries)

    return new_Boundaries,dists

        
def getWidth(Boundary,*args):
    """
    Given the Boundary points of a set D at each location x* in xRange, get the width of the set at the Intersection of (x*,y) and D, it
    is assumed that this intersection is a Finite Union of Intervals.
    """
    Width =[]
    for b in Boundary:
        Width.append(SumBound(b))
        
    return Width

def SumBound(B):
    """
    B is a list of tuples [(a1,a2),...], that may contain null elements specified by (0,), get the sum of the difference between the elements of each tuple.
    """
    Sum = 0
    for Tuple in B:
        if len(Tuple)==1:
            break
        else:
            Sum+=(Tuple[1] - Tuple[0])
    return MyFloat(Sum)


def Convert(X,Y):
    """Fuse the elements of the sublists of X and Y in a tuple and append it to the array 'new_array'
        record the number of elements in each sublist in X(Y) in the list dist
        @param X   a List of lists contanining floating point elements
        @param Y  a List of lists contanining floating point elements
    """
    new_array = []
    dist = []
    for subarray in range(len(X)):
        dist1 = len(X[subarray])
        for i in range(dist1):
            new_array.append(MyTuple([MyFloat(X[subarray][i]),MyFloat(Y[subarray][i])]))
        dist.append(dist1)
    if len(dist)>1:
        dist = MyTuple(dist)
    return new_array,dist

def ConstructArray(List_Arrays):
    """Reads in a list of arrays, calculates the one with the largest number of elements an storage it in the value Max_length
       completes all the other arrays to that number by filling the missing elements with NaN, after that appends all of them
       in the list Array_handler, convert it to a numpy array and transpose it.
       @param List_Arrays :  a List of lists 
       """
    Max_Length = max([len(Array) for Array in List_Arrays])
    Array_handler = []
    for array in List_Arrays:
        length = len(array)
        n = Max_Length - length
        for i in range(n):
            array.append('NaN')
        array.append('NaN')
        Array_handler.append(array)
    return np.array(Array_handler).transpose()

def FormatEq(zeq):
    return MyTuple([MyFloat(x) for x in zeq])

def FormatStab(zstab):
    Eigenvals = []
    for eigenval in zstab:
        x = eigenval.real
        y = eigenval.imag
        Eigenvals.append(MyInnerTuple([MyFloat(x),MyFloat(y)]))
    return MyTuple(Eigenvals)
        
def inList(yPoint,ListofTuples):
    for T in ListofTuples:
        if len(T)==1:
            return False
        Ylog= np.log10(yPoint)
        if T[0]<Ylog and Ylog<T[1]:
            return True
    return False

def getRegion(xIndex,yPoint,PPos1,CPos1,PPos2,CPos2,DPos,fDict,args):
    r"""
    returns the zone in which a given yPoint belongs, the zones are delimited by the following criterions:

    .. math::
    

       &1 \iff  y \in Z(I_{C4}) \cap Z(I_{C5})\\
       &2 \iff  (y \in Z(I_{C4}) \cap Z(I_{C3}) \land y \notin Z(I_{C5})) \lor ( y \in Z(I_{C3}) \setminus Z(I_{C5})  \land y \notin Z(I_{C2})) \lor (y \in Z(I_{C3}) \setminus Z(I_{C5}) \cap Z(I_{C2}) \setminus Z(I_{C4}) \land (I_3(y) > I_2(y)) ) \\
       &3 \iff (y \in Z(I_{C5}) \cap Z(I_{C2}) \land y \notin Z(I_{P4})) \lor (y \in Z(I_{C2}) \setminus Z(I_{C4}) \land y \notin Z(I_{C3}) ) \lor (y \in Z(I_{C2})\setminus Z(I_{C4}) \cap Z(I_{C3})\setminus Z(I_{C4})) \land (I_2(y) > I_3(y) ) \\
       &4 \iff \mbox{any of the previous criterions are false}\\

 
       
    """
    tP1 = inList(yPoint,PPos1[xIndex])
    tC1 = inList(yPoint,CPos1[xIndex])
    tP2 = inList(yPoint,PPos2[xIndex])
    tC2 = inList(yPoint,CPos2[xIndex])

    
    return setr1(tP1,tP2,tC1,tC2,fDict,args)

def setr1(tP1,tP2,tC1,tC2,fDict,args):
    """
    from a set of boolean values returns an integer which is the index of a particular region
    """

    if tP2 and tC2:
        return 1
    elif tP2 and tC1 and tP1:
        return 2
    elif tC2 and tP1 and tC1:
        return 3
    elif tP1 and not tC1:
        return 2
    elif tC1 and not tP1:
        return 3
    elif tC1 and tP1:
        a = fDict['I_C_s2'](*args)
        b = fDict['I_P_s3'](*args)
        if a > b :
            return 3
        else:
            return 2
    else:
        return 4

                
def getRegion2(xIndex,yPoint,PPos1,CPos1,PPos2,CPos2,DPos,fDict,args):
    r"""
    Same ase the getRegion functions, except that account for the cases where exists a combinations of (x,y) such that D<0; if that is the case
    the criterions change to:

    .. math::
     
       1 \iff y \notin ( Z(I_{C4}) \lor Z(I_{C5})) \\
       2 \iff ((y \in (Z(I_{C4}) \cap Z(I_{C5}))) \land (I_5(y) > I_4(y))) \cup (2-getRegion1) \\
       3 \iff ((y \in (Z(I_{C4)) \cap Z(I_{C4}))) \land ( I_4(y) > I_5(y))) \cup (3-getRegion1) \\
       4 \iff \mbox{if none of the above holds}
    

    In words this assumes that the noninvasible equilibrium state is picked as the final state or if there is no such state then it is picked the biggest one that is \emph{less likely} to be invaded, that is the one for which the invasibility function of the invader gets a smaller value. In the case for which there is more than one noninvasible equilibrium state it is chosen the one that is more likely to arise, that is the one for which the invader invades \emph{faster}
    """
    tP1 = inList(yPoint,PPos1[xIndex])
    tC1 = inList(yPoint,CPos1[xIndex])
    tP2 = inList(yPoint,PPos2[xIndex])
    tC2 = inList(yPoint,CPos2[xIndex])
    tD = inList(yPoint,DPos[xIndex])

    if tD:
        return setr1(tP1,tP2,tC1,tC2,fDict,args)
    else:
        if not(tP2) and not(tC2):
            return 1
        elif tP2 and tC2 and tP1 and tC1:
            a = fDict['I_C_s5'](*args)
            b = fDict['I_P_s4'](*args)
            if a > b:
                return 3
            else:
                return 2
        elif tP2 and tC1 and tP1:
            return 2
        elif tC2 and tP1 and tC1:
            return 3
        elif tP1 and not tC1:
            return 2
        elif tC1 and not tP1:
            return 3
        elif tC1 and tP1:
            a = fDict['I_C_s2'](*args)
            b = fDict['I_P_s3'](*args)
            if a > b:
                return 3
            else:
                return 2
        else:
            return 4
            
        
def GetUnion(Bound1,Bound2):
    r"""
    get the boundary points of the union of two sets A and B, having as input the set of boundary points at each :math:`x \in xRange`
    """
    Union = []
    for x in range(len(Bound1)):
        U = FindUnion(Bound1[x],Bound2[x])
        Union.append(U)

def FindUnion(BP1,BP2):
    if len(BP1[0])==1 and len(BP2[0]) == 1:
        return [(0,)]
    else:
        Int = FindIntersection(BP1,BP2)
        CompInt = GetComplement(Int,Min,Max)
        
        A1 = FindIntersection(CompInt,BP1)
        A2 = FindIntersection(CompInt,BP2)
        U = FormatUnionSet(sorted(A1+A2+Int))
        return U


def FormatUnionSet(L):
    """
    Glue together the intervals of the sorted List L,sorted in increasing order, if its limits coincide . Each interval is represented by a Tuple . 
    It uses a recurse algorithm in which we start at the leftmost position L[0], proceed to the right gluing together all the intervals to which it coincide ,i.e. L[0][1] = L[i][0], it keeps updating L[0] until it don't found more, and repeat the same procees starting from that position. It terminates when there are no more intervals to add.
    """
    newL = []
    n = range(L)
    while i < n :
        if i == n-1:
            newL.append(L[i])
        else:
            k = 1
            while i+k < n and  rangeL[i+k][0] == L[i][1]:
                L[i][1] = L[i+k][1]
                k+=1

            
            newL.append(L[i])
            i += k

    return newL
    

    
    
def GetIntersection(Bound1,Bound2):
    r"""
    Get the boundary of the intersection of two sets A and B, by having as input the set of boundary points at each :math:`x \in xRange`
    """
    Intersection=[]
    for x in range(len(Bound1)):
        if len(Bound1[x][0])>1:
            Int = FindIntersection(Bound1[x],Bound2[x])
            Intersection.append(Int)
        else:
            Intersection.append([(0,)])
    return Intersection

def FindIntersection(BoundaryPoints1,BoundaryPoints2):
    r"""
    find the intersection of two sets of intervals A and B whose limit points are stored in BoundaryPoints1 and BoundaryPoints2 respectively.
    It is based on the fact that :

    .. math::

       (\cup_{I \in A} I) \cap (\cup_{J \in B} J) = \cup_{J \in B} (\cup_{I \in A} (J \cap I))


    """
    IntPoints=[]
    if len(BoundaryPoints2[0])>1: 
        for i in range(len(BoundaryPoints2)):             
            NewBPoints= findIntPoints(BoundaryPoints1,BoundaryPoints2[i])
            addPoints(IntPoints,NewBPoints)     
    else:
        IntPoints=[(0,)]
        
    IntPoints = FormatStart(IntPoints)
    
    return IntPoints

def findIntPoints(BoundaryPoints1,BoundaryPoint2):
    r"""
    returns the boundary points of :math:`(\cup_{I \in A} (J \cap I))` for a particular interval :math:`J` whose boundary points are given by BoundaryPoint2 and the boundary points of the intervals I are given by the list BoundaryPoints1
    """
    Points=[]
    for bp in BoundaryPoints1:
        Int = Intersection(BoundaryPoint2,bp)
        addPoints(Points,Int)
    
    #Format Points
    Points = FormatStart(Points)
      
    return Points

def FormatStart(Points):
    """
    Format the list of Points to account that [(0,)] elements could have been added to by the first iteration, that is the first interval(s) could have a null intersection. But the overall intersection is not, in that cases it deletes the first element.

    """
    if len(Points) > 1:
        if len(Points[0]) == 1:
            Points = Points[1:]
    return Points
        
        
    
def Intersection(bp1,bp2):
    r"""
    For intervals I and J whose boundary points are given by bp1 and bp2 respectively , returns the boundary points for :math:`I \cap J`
    """
    if bp1[0]<bp2[0]:
        return IntervalInt(bp1,bp2)       
    else:
        return IntervalInt(bp2,bp1)
    
def IntervalInt(bpL,bpH):
    r"""
    Find the boundary points of :math:`(a,b) \cap (c,d)` for which :math:`a \leq c`
    """
    if bpH[0]<bpL[1]:
        if bpH[1]<bpL[1]:
            return [bpH]
        else:
            return [(bpH[0],bpL[1])]
                
    else:
        return [(0,)]
        
def addPoints(IntPoints,NewBPoints):
    """
    Add NewBPoints **N** to the List IntPoints **I**, only add :math:`x \in N` if :math:`|x[0]|>1 \lor |\mathbf{I}| = 0`
    """
    if len(NewBPoints[0])==1:
        if len(IntPoints)>0:
            pass
        else:
            IntPoints.append(NewBPoints[0])
    else:
        for i in range(len(NewBPoints)):
            IntPoints.append(NewBPoints[i])
## Output Clases

def constructFooter(params):
    s = 'params used'
    for param in params:
        s += ' ' + param + '='+str(params[param])
    s= s[:]
    return s

import csv
class Data(object):
    """Class that stores a dataset in the format of a n x n matrix """
    def __init__(self,data):
        self.data = data 
     
    def setData(self,Data):
        self.data = Data
        
    def ncols(self):
        return len(self.data[0])
    
    def nrows(self):
        return len(self.data)
    def getrow(self,index):
        return self.data[index]
    def TransformtoFloats(self):
        #Transform dataset entries to floats
        for i in range(self.nrows()):
            for j in range(len(self.getrow(i))):
                if self.data[i,j] != 'NaN':
                    self.data[i,j] = [ float(item) for item in self.data[i,j][1:-1].split(':')]
            
    def reshape(self):
        self.data.reshape((self.nrows(),self.ncols()))



class OutputInvData(Data):
    """ Creates an abstraction of a csv table, with a header and footer.
        data stores all the body of the table
        distribution is a tuple whose length referes to the number of parts which has been put together
        in the column, and each component of it referees to the number of items in each part"""
    def __init__(self,data,header,footer,xFocusAndSep,distribution):
        self.data = data
        self.header = header
        self.footer = footer
        self.distribution = distribution
        self.xFocusAndSep = xFocusAndSep

    def setHeader(self,header):
        self.header = header
    
    def setFooter(self,footer):
        self.footer = footer
                
    
    def WriteWidths(self,direction,delimiter):
        with open(direction,'wb') as csvfile:
            writer = csv.writer(csvfile,delimiter=delimiter,quotechar='|',quoting= csv.QUOTE_MINIMAL)
            writer.writerow(self.header)
            writer.writerows(self.data)
            writer.writerow(self.footer)
            writer.writerow(self.xFocusAndSep)

    def WriteInvasibility(self,direction,delimiter):
        """Make use of the csv module to produce an csv file with the data contained in self.data """
        with open(direction,'wb') as csvfile:
            writer = csv.writer(csvfile,delimiter=delimiter,quotechar='|',quoting= csv.QUOTE_MINIMAL)
            writer.writerow(self.header)
            writer.writerows(self.data)
            writer.writerow(self.distribution)
            writer.writerow(self.footer)
            writer.writerow(self.xFocusAndSep)
            
    def WriteEquilibrium(self,direction,delimiter):
        with open(direction,'wb') as csvfile:
            writer = csv.writer(csvfile,delimiter=delimiter,quotechar='|',quoting= csv.QUOTE_MINIMAL)
            writer.writerows(self.data)
            writer.writerow(self.footer)
            writer.writerow(self.xFocusAndSep)












##Dynamics class
class Dynamics(BSR):
    def __init__(self,workingData,initCondition,finalTime,separation,K_CP,K_RC,m_P):
        BSR.__init__(self,workingData.getParams(),workingData.getmode(),workingData.getmassLims())
        self.initCondition = initCondition
        self.finalTime = finalTime
        self.separation = separation
        self.K_CP = K_CP
        self.K_RC = K_RC
        self.m_P = m_P
        self.dR = ''
        self.dC = ''
        self.dP = ''
        self.fDict = workingData.getfDict()
        self.Kpoints = 1e4
    def setParamVals(self,K_CP,K_RC,m_P):
        """ Set the values for the three key paramaters of the model , the size ratios and the body mass """
        self.K_CP = K_CP
        self.K_RC = K_RC
        self.m_P = m_P
            
    def setinitCondition(self,initCondition):
        """ Specifices the given initial condition from which to start the simulation """
        self.initCondition = initCondition
        
    def makeinitCondition(self,case):
        
        """ se the init condition , depending on the scenario , in the Invasibility by P to C-R the initial condition
        is the equilibrium of the latter two in isolation , a similar situation is for the invasibility of C to P-R(labeled 
        scenario 5)"""
        if self.mode == "RM":
            if case == "Inv P4":
                self.initCondition[0] = self.fDict['R_eq_s2RM'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['C_eq_s2RM'](self.K_CP,self.K_RC,self.m_P)
            elif case == "Inv P5":
                self.initCondition[0] = self.fDict['R_eq_s3RM'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['P_eq_s3RM'](self.K_CP,self.K_RC,self.m_P)
        elif self.mode == "LV":
            if case == "Inv P4":
                self.initCondition[0] = self.fDict['R_eq_s2'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['C_eq_s2'](self.K_CP,self.K_RC,self.m_P)
            elif case == "Inv P5":
                self.initCondition[0] = self.fDict['R_eq_s3'](self.K_CP,self.K_RC,self.m_P)
                self.initCondition[1] = self.fDict['P_eq_s3'](self.K_CP,self.K_RC,self.m_P)
                               
            
    def setfinalTime(self,finalTime):
        """ set the time until when to stop the simulation """
        self.finalTime = finalTime
    def setseparation(self,separation):
        self.separation = separation
            
    def setDynamicFunction(self):
        """ input the corresponding values for the paramters K_CP, K_RC and m_P and convert the dynamical functions dR, dC 
        and dP in three-argument functions, just depending on the value of the biomass densities R, C and P"""
             
        args = np.array([self.R,self.C,self.P,self.K_RC,self.K_CP,self.m_P])
        if self.mode == "RM":
            self.dR = lambdify((self.R,self.C,self.P),self.fDict['dRRM'](*args))
            self.dC = lambdify((self.R,self.C,self.P),self.fDict['dCRM'](*args))
            self.dP = lambdify((self.R,self.C,self.P),self.fDict['dPRM'](*args))
        else:
            self.dR = lambdify((self.R,self.C,self.P),self.fDict['dRLV'](*args))
            self.dC = lambdify((self.R,self.C,self.P),self.fDict['dCLV'](*args))
            self.dP = lambdify((self.R,self.C,self.P),self.fDict['dPLV'](*args))
            
       
    def DynamicFunction(self,X,t):
        args = np.array([X[0],X[1],X[2]])
        
        return np.array([self.dR(*args),self.dC(*args),self.dP(*args)])
    
    
        
    def Simulate(self):
        """ simulation routine using the Odeint solver from the scipy optimize package , Odeint is a python implementation
        of the ODEPACK package in FORTRAN which uses a multi-step solver in the non stiff case """
        self.setDynamicFunction()
        t = np.arange(0,self.finalTime,self.separation)
        return odeint(self.DynamicFunction,self.initCondition,t,())
    
    
    def runSimulationSimK(self,case,massLims,lowKLims,upKLims,initDirection):
        for massIndex in xrange(len(massLims)):
            
            Krange = 10**(np.linspace(lowKlims[massIndex],upKLims[massIndex],self.Kpoints))
            Kdata = []
            for K in Krange:
                self.setParamVals(K,K,massLims[massIndex])
                self.setinitCondition(case)
                Run = self.Simulate()
                Kdata.append(Run)
            WriteData(Kdata,initDirection+str(massIndex)+".csv") 


