from Baseclass import *

class InvBoundaries(BSR):
    """ 
    Class that stores and computes the values for the invasibility boundaries of the distinct scenarios, it receives as initial
    input the values of the parameters and the invasibility functions 
    """
    def __init__(self,workingData):
        
        BSR.__init__(self,workingData.getParams(),workingData.getmode(),workingData.getxLims(),workingData.ksim)
        self.xRange = workingData.setxRange()
        self.yRange = workingData.setyRange()
        self.UpGuess =  8
        self.LowGuess = -13.005
        self.guessSep = 0.01
        self.InvFunctions=workingData.setInvFunctions()
        self.InvFunctionsNames = ['I_C_s2','I_P_s3','I_P_s4','I_C_s5']
        self.InvBounds={}
        self.UnEditedInvBounds={}
        self.fDict = workingData.getfDict()
        self.distance = 1e-12
        self.PositiveBoundaries={}
        self.currentMass = self.params['massP']
        self.xFocus = workingData.xFocus
        self.xFocus_sep = workingData.xFocus_sep
        self.DBound = 0
        self.Intersections={}
        self.Widths={}
        self.SZ={}
        self.ZBounds={}
        self.Footer = constructFooter(self.params)
        self.EfDif = self.params['e1']*self.params['e3'] - self.params['e2']
            
    def setUpGuess(self,Guess):
        self.UpGuess = Guess
            
    def UpdateMass(newMass):
        self.currentMass = newMass
        
    def setAndWriteInvBoundaries(self,Header,Direction,Direction2):
        """
        * Find the zero boundaries for each of the functions specified in self.InvFunctions
        * Write them to a csv file specified in Direction.
        
        Header refers to the first row of the csv file.
        """

        self.setInvBoundaries()
        self.writeInvasibilityValues(Header,Direction)
        self.writeUneditedInvValues(Header,Direction2)
        
        
    def setInvBoundaries(self):
        """ 
        For each of the functions present in the InvFunctions list, computes the invasibility boundaries by means of the
        Scipy.Brentq numerical method.
        For each of the invasibility functions it returns a dict object with x and y keys indicating the zeros(size ratio values delimiting
        the boundaries) for the function at each x (mass), the values of x and y are a list of lists. the total number of sublists denotes the maximum
        number of zeros found at any location x.
        """
        searchRange = 10**(np.arange(self.LowGuess,self.UpGuess,self.guessSep))
        xRange = self.xRange
        if self.EfDif<0:
            self.InvFunctionsNames+=['D']
            self.InvFunctions['D'] = self.fDict['D']


        for invfunc in self.InvFunctionsNames:
            bound_dict = self.InvBoundary(invfunc,searchRange,xRange)
            self.InvBounds[invfunc] = bound_dict           
    

    def InvBoundary(self,Invfunc,searchRange,xRange):
        """Find the inv boundaries(zeros of the invasibility function) using the brentq method from the SciPy package,
        input arguments, the invasibility function, the limits for the interval to look for zeros(depending on the x values) """
        
        f = self.InvFunctions[Invfunc]
        K = Get_bounds2(f,Get_roots,xRange,searchRange,additionalPar= self.currentMass,k_sim = self.ksim)
        x,y= procce_(K,xRange)
        self.UnEditedInvBounds[Invfunc]=K
        return {'x':x,'y':y}
    
    def getBounds(self):
        """ return the dictionary containing all the invasibility boundaries"""
        return self.InvBounds
    

    def writeUneditedInvValues(self,Header,direction,delimiter=','):
        if self.EfDif<0:
            H = Header + ['D']
        else:
            H = Header
        
        newData = EditData(self.InvFunctionsNames,self.UnEditedInvBounds)
        newData = ConstructArray(newData)
        OutputFile = OutputInvData(newData,H,[self.Footer],[MyTuple(self.xFocus),self.xFocus_sep],'')
        OutputFile.WriteWidths(direction,delimiter)
    def writeInvasibilityValues(self,Header,direction,delimiter=','):
        """Write data into a csv file specified in direction
        @param Boundaries a list of 2-elements lists which each of them stores the X and Y coordinates
        of the invasibility boundaries computed in the analysis
        @param the list of params used to compute the boundaries which will be used in the footer of the csv file
        @direction the system direction where the file is going to be stored
        @Header the first row of the csv file """
       

        if self.EfDif<0:
            H = Header + ['D']
        else:
            H = Header
        new_Boundaries,dists = FormatZones(self.InvFunctionsNames,self.InvBounds)
        Footer = self.Footer
        OutputFile = OutputInvData(new_Boundaries,H,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],dists)
        OutputFile.WriteInvasibility(direction,delimiter)

       
    def setPositiveBoundaries(self):
        """ Creates a dict storying the boundary points of the set in which each of the criterions is satisfied """
        for Inv in self.UnEditedInvBounds.keys():
            M,P = RME.GetGuessBounds(self.UnEditedInvBounds[Inv],self.xRange,self.fDict[Inv],self.UpGuess,self.LowGuess,self.ksim,self.currentMass)
            
            self.PositiveBoundaries[Inv] = P
            
    def WriteWidths(self,Direction):
        """
        Write the Width dictionary to a file specified in Direction.
        """
        Header = self.Widths.keys()
        data = FormatWidths(Header,self.Widths)
        Footer = self.Footer
        OutputFile = OutputInvData(data,Header,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],"")
        OutputFile.WriteWidths(Direction,',')
        

    def setAndWriteWidthsZones(self,DirectionW,DirectionZB):
        """
        * For each x in xRange and any Invasibility function f calculates the sum of the length of all the intervals which are in the Positive region of f
        * Get the Positive Boundaries for all the target zones in the analysis
        * Write both of the above results to a file whose direction is specified in the arguments DirectionW and DirectionZ.
        """
        self.setWidthsAndZones()
        self.WriteWidths(DirectionW)
        self.WriteZonesBounds(DirectionZB)
        #self.WriteZones(DirectionZones)
    
    def setWidthsAndZones(self):
        """
        * Calculates the Boundaries of the Zones described in the Study, which are intersection of the Positive Regions of a subset of the Invasibility functions
        * Calculates the widths of each of the Zones
        * Convert them to the format used for the Invasibility functions and save them in the ZBounds dict
        """

        self.setIntersections()
        self.setZonesBoundaries()
        #self.setZones()
        self.setWidths()

    def setIntersections(self):
        r"""
        Find the boundary points of the Intersections of the Positive Regions of some invasibility functions, which delineate the zones expressed in the analysis.
        For example :math:`Z(I_{C4}) := Z(I_{C2}) \cap PR_{4}` where :math:`PR_{4}` is the set for which the Invasibility function I_P_s4 is positive.
        We assume that one dimensional subset of the set are always a union of \emph{connected} sets.
        """
        
        I1 = GetIntersection(self.PositiveBoundaries['I_C_s5'],self.PositiveBoundaries['I_P_s3'])
        I2 = GetIntersection(self.PositiveBoundaries['I_C_s2'],self.PositiveBoundaries['I_P_s4'])
        I3 = GetIntersection(I1,I2)
        I4 = GetIntersection(self.PositiveBoundaries['I_C_s5'],self.PositiveBoundaries['I_P_s4'])


        self.Intersections['Z(IC5)'] = I1
        self.Intersections['Z(IC4)'] = I2
        self.Intersections['MutualInv'] = I3

        
        if self.EfDif < 0 :

            I5 = GetIntersection(GetComplement(self.PositiveBoundaries['I_C_s5'],self.LowGuess,self.UpGuess),GetComplement(self.PositiveBoundaries['I_P_s4'],self.LowGuess,self.UpGuess,))
            I6 = GetIntersection(I5,GetComplement(self.PositiveBoundaries['D'],self.LowGuess,self.UpGuess,))
            I7 = GetIntersection(I4,self.PositiveBoundaries['D'])

            self.Intersections['UnstableCoexistence'] = I6
            self.Intersections['PotentialStableCoexistence'] = I7
        else:
            self.Intersections['PotentialStableCoexistence'] = I4

        
        self.findStabilityBounds()

    def findStabilityBounds(self):
        f = self.fDict['hd2']
        BSearchRange = self.Intersections['PotentialStableCoexistence']

        K = Get_boundsASR(f,Get_roots,self.xRange,BSearchRange,additionalPar=self.currentMass,k_sim = self.ksim)
        self.Intersections['StableCoexistence'] = RME.GetGuessBoundsGC(K,self.xRange,f,BSearchRange,self.ksim,self.currentMass)[1]
    

    def setWidths(self):
        """
        Find the widths of each of the Zones 
        """
        F(self.Widths,getWidth,self.PositiveBoundaries,self.Intersections,self.EfDif)
    def setZonesBoundaries(self):
        """
        Find the boundary of the Zones and format them to the same data structure used in the computation of the Invasibility boundaries
        """
        F(self.ZBounds,GetPositiveBoundaries,self.PositiveBoundaries,self.Intersections,self.EfDif,self.xRange)

        
    def findZones(self):
        """
        From each of the Boundaries, create a 2D array using xRange and yRange and return a 1 0  array , each 1 located at a position in which
        the point is interior to the Positive Region delimited by the Boundary
        """
        F(self.Zones,GetPositiveRegions,self.PositiveBoundaries,self.Intersection,self.xRange,self.yRange)


    def WriteZonesBounds(self,Direction):
        """
        Write the SZBounds into a file whose pointer is specified in Direction.
        """
        Header = self.ZBounds.keys()
        data,dist = FormatZones(Header,self.ZBounds)
        Footer = self.Footer
        OutputFile = OutputInvData(data,Header,[Footer],[MyTuple(self.xFocus),self.xFocus_sep],dist)
        OutputFile.WriteInvasibility(Direction,',')


        

        

