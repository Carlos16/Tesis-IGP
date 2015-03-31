from AuxiliarClases import *
from bounds import *
import RMEquibrium as RME

########## Auxiliary functions######################        
def GetPositiveRegions(Bounds,xRange,yRange):
    X,Y = np.meshgrid(xRange,yRange)
    m,n = X.shape
    Z = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            Z[i,j] = int(inList(Y[i,j],Bounds[j]))
    return Z

            
def F(D,F1,D1,D2,*args):
    D['Z(IC2)'] = F1(D1['I_C_s2'],*args)
    D['Z(IC3)'] = F1(D1['I_P_s3'],*args)
    D['Z(IC4)'] = F1(D2['Z(IC4)'],*args)
    D['Z(IC5)'] = F1(D2['Z(IC5)'],*args)
    D['MutualInv'] = F1(D2['MutualInv'],*args)    

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
    else:
        if len(dist) == 0:
            dist = MyTuple((0,0))
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
    return np.array(Array_handler,dtype = object).transpose()

        
def inList(yPoint,ListofTuples):
    r"""
    Returns True if :math:`y \in T , T \in LoT` , where *T* is a tuple and *LoT* is a list of tuples
    """
    for T in ListofTuples:
        if len(T)==1:
            return False
        Ylog= np.log10(yPoint)
        if T[0]<Ylog and Ylog<T[1]:
            return True
    return False

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
