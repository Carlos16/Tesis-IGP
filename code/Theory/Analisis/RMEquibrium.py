from bounds import *



def FindCoexistenceRegion(self):
    self.EquilibriumBoundaries['R'] = self.FindRpositiveRegion()
    self.EquilibriumBoundaries['C'] = self.FindCpositiveRegion()
    self.EquilibriumBoundaries['P'] = self.FindPpositiveRegion()


def FindRPositiveRegion(self):
    UpGuessRange,LowGuessRange,massRange  = self.MakeRangeGuess(False,'Discriminant')
    DisBounds = Get_bounds(self.fDict['Discriminant'],Get_roots,massRange,UpGuessRange,LowGuessRange,self.guessSep)
    R1Bounds = [[massRange,Get_bounds(self.fDict['R1'],Get_roots,massRange,UpGuessRange,LowGuessRange,self.guessSep)]]
    R2UpGuessRange,R2LowGuessRange,R2massRange = GetGuessBounds(DisBounds,massRange,self.fDict['Discriminant'],self.UpGuess,self.LowGuess,self.distance)
    R2Bounds =[]
    R3Bounds =[]
    for j in range(len(R2UpGuessRange)):
        R2Bounds.append([massRange[j],Get_Bounds(self.fDict['R2'],Get_roots,R2massRange[j],R2UpGuessRange[j],R2LowGuessRange[j],self.guessSep)])
        R3Bounds.append([massRange[j],Get_Bounds(self.fDict['R3'],Get_roots,R2massRange[j],R2UpGuessRange[j],R2LowGuessRange[j],self.guessSep)])
    return {'R1':R1Bounds,'R2':R2Bounds,'R3':R3Bounds}



def FindpositiveRegion(self,Pred):
    Bounds={}
    REquilibrium = self.EquilibriumBoundaries['R']
    for Rpos in REquilibrium.keys():
        Boundaries = []
        for Bound in REquilibrium[Rpos]:
                UpGuessRange,LowGuessRange,massRange =GetGuessBounds(Bound[1],Bound[0],self.fDict[Rpos],self.UpGuess,self.LowGuess,self.distance)
                for j in range(len(UpGuessRange)):
                    Boundaries.append([massRange[j],Get_Bounds(self.fDict[Pred][Rpos],Get_roots,massRange[j],UpGuessRange[j],LowGuessRange[j],self.guessSep)])
        Bounds[Rpos]=Boundaries
    return Bounds




#Delimit Positive Boundaries and Return Lists

def ExtractPoints(xRanges,BoundaryPoints):
    xGuess=[[]]
    UpGuess=[[]]
    LowGuess=[[]]
    for x in range(len(xRanges)):
        if len(BoundaryPoints[x][0])>1:
           for p in range(len(BoundaryPoints[x])):
               try:
                   xGuess[p].append(xRanges[x])
                   UpGuess[p].append(BoundaryPoints[x][p][1])
                   LowGuess[p].append(BoundaryPoints[x][p][0])
               except:
                
                   xGuess+=[[xRanges[x]]]
                   UpGuess+=[[BoundaryPoints[x][p][1]]]
                   LowGuess+=[[BoundaryPoints[x][p][0]]]
    return UpGuess,LowGuess,xGuess
                
            

        
def toDict(xRanges,BPoints):
    BPointDict={}
    for x in range(len(xRanges)):
        BPointDict[xRanges[x]] = BPoints[x]

    return BPointDict


def GetGuessBounds(Arr,Arr2,f,defaultUpGuess,defaultLowGuess,ksim,arg1):
    """ Code for delimiting the boundaries of the positive region of a given function"""
    xRanges=[]
    BoundaryPoints=[]
    for i in range(len(Arr)):
        AddToBoundaryPoints(Arr[i],Arr2[i],xRanges,BoundaryPoints,f,defaultUpGuess,defaultLowGuess,ksim,arg1)

    return xRanges,BoundaryPoints


def GetGuessBoundsGC(Arr,Arr2,f,DomainBounds,ksim,arg1):
    xRanges=[]
    BoundaryPoints=[]
    for i in range(len(Arr)):
        AddtoBoundaryPointsG(Arr[i],Arr2[i],xRanges,BoundaryPoints,f,DomainBounds[i],ksim,arg1)
    return xRanges,BoundaryPoints

def AddtoBoundaryPointsG(BreakPoints,xPoint,xList,yList,f,DomainBoundP,ksim,arg1):
    wList =[]
    xList.append(xPoint)
    if len(DomainBoundP[0])>1:
        for i in range(len(DomainBoundP)):
            uList=[]
            AddToBoundaryPoints(BreakPoints[i],xPoint,[],uList,f,DomainBoundP[i][1],DomainBoundP[i][0],ksim,arg1,eps=1e-10,eps2=1e-10)
            if len(uList[0][0])>1:
                wList+=uList[0]

        if len(wList)==0 :
            yList.append([(0,)])
        else:
            yList.append(wList)
    else:
        yList.append([(0,)])
                            

def AddToBoundaryPoints(BreakPoints,xPoint,xList,yList,f,maxY,minY,ksim,arg1,eps=1e-03,eps2=0.):    
    n = len(BreakPoints)
    xList.append(xPoint)
    if n>0:
        if isNegative(f,xPoint,minY+eps2,ksim,arg1):
            AddPoints(BreakPoints,yList,n,maxY,minY,eps,Type=-1)
        else:
            AddPoints(BreakPoints,yList,n,maxY,minY,eps,Type=+1)
    if n==0:
        if not isNegative(f,xPoint,minY+eps2,ksim,arg1):
            yList.append([(minY-eps,maxY+eps)])
        else:
            yList.append([(0,)])


def AddPoints(BreakPoints,yList,n,maxY,minY,eps,Type):
    """ From a list of zeros determine the ones that delimit a pass from positive to negative and store them in a list, i.e. find the boundary points
    of the positive region"""
    initPos = np.sign(1 - Type)
    BreakPoints = [minY-eps]+BreakPoints+[maxY+eps]
    Points=[]
    for i in range(initPos,n+1,2):
        Points.append((BreakPoints[i],BreakPoints[i+1]))
    yList.append(Points)



    
        

def isNegative(f,xPoint,yPoint,ksim,arg1):
    if ksim:
        return f(10**yPoint,xPoint) <0
    else:
        return f(10**yPoint,xPoint,arg1)<0


    

        
        
    

        
            
            
        
        
        
    
    

    
    
    
        
