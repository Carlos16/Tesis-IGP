import numpy as np

class Point(object):
    """ Class that stores the information of a two dimensional point in a semilog scale , y is in log scale"""
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def dist(self,Point2):
        """Computes the distance between two points in the actual scale"""
        return np.sqrt((self.x-Point2.x)**2 + (Point2.y-self.y)**2)
    
    def logdist(self,Point2):
        """Transform x to log scale and computes the euclidean distance """
        return np.sqrt((np.log10(self.x)-np.log10(Point2.x))**2 + (Point2.y-self.y)**2)

    def getx(self):
        return self.x
        
    def gety(self):
        return self.y

    def getCoordinates(self):
        return (x,y)
    
    def __str__(self):
        return '('+str(self.x) + ',' + str(self.y) +')'
    
    def __repr__(self):
        return self.__str__()

    def __eq__(self,P):
        if P.x == self.x and P.y ==self.y :
            return True
        else:
            return False
        
    def __lt__(self,P):
        if P.x >self.x:
            return True
        else:
            return False
    
class Path(object):
    """ Class that stores the information of a two dimensional set of points , each point is a tupple, it divides a path in its coneccted components, which are
    divided based on a threshold distance that consecutive points(in term of the x axis) must satisfy, each division is called a subpath and the ends of each of t
    the subpaths are stored in the identifiers variable.
    It also provides methods for creating a path,fusing two paths and extracting its points 
    for creation it recieves a stream of points for the X and Y components"""
    def __init__(self,subPathsX,subPathsY):
        self.subpaths = self.CreatePath(subPathsX,subPathsY)
        self.Identifiers =[[subpath[0],subpath[-1]] for subpath in self.subpaths]
        self.Dist= 0.2
        self.xDist = 0.025
        
    def setDist(self,D):
        self.Dist = D
    def getSubpath(self,i):
        """output the subpath at index i"""
        return self.subpaths[i]
    def UpdateSubpath(self,path,i):
        """Replaces the subpath at index i with a new subpath path"""
        self.subpaths[i] = path
    def deleteSubpath(self,i):
        """Delete the subpath stored in the index i"""
        self.subpaths.pop(i)
        self.Identifiers.pop(i)
        
    def getSubPathXY(self,index):
        P = self.subpaths[index]
        subX = [p.x for p in P]
        subY = [p.y for p in P]
        return subX,subY
    def CreatePath(self,subPathsX,subPathsY):
        """Create a list of sublists of points, in which each of the elements is a point.
        It recieves two list of lists of scalars which represents x and y coordinates respectively as input."""
        Path = []
        for n in range(len(subPathsX)):
            A = [Point(subPathsX[n][i],subPathsY[n][i]) for i in xrange(len(subPathsX[n]))]
            Path.append(A)
        return Path
    
    def AddPath(self,Path):
        """Combine two Paths, it adds a second Path X to the existing Path Y, it adds each of the subpaths of X one at the time"""
        for index in range(len(Path.subpaths)):
            self.addSubpath(Path,index)

        
    def FormatPoints(self):
        Points_indexes = []
        
        for i in range(len(self.subpaths)):
            if len(self.subpaths[i]) ==1:
                Points_indexes.append(i)
        Point_list = []
        Point_ident = []
        for i in Points_indexes:
            Point_list.append(self.subpaths[i][0])
            Point_ident.append(self.Identifiers[i][0])
            #self.deleteSubpath(i)
        
        if len(Points_indexes)>0:
            self.subpaths[Points_indexes[0]-1]+=Point_list
            self.Identifiers[Points_indexes[0]-1][1] = Point_ident[-1]
        self.RemovePoints()
    
    def RemovePoints(self):
        SubPathHolder = []
        IdentHolder = []
        for i in range(len(self.subpaths)):
            if len(self.subpaths[i])>1:
                SubPathHolder.append(self.subpaths[i])
                
                IdentHolder.append(self.Identifiers[i])
        self.subpaths = SubPathHolder
        self.Identifiers = IdentHolder
    
    def addSubpath(self,Path,i):
        """Add the subpath of index i,
          It uses the identifiers of the subpath to find the closest subpaths of self.Path,if they are belong they are belong the threshold distance, and fuse each of the extreme elements of the subpath to the respective subpaths(ends). if none of the subpaths of self.Path are below the threshold it just add the subpath as a new subpath of self.Path   

        """
        Id = Path.Identifiers[i]
        p,r = Match(Id,self.Identifiers,self.Dist)    #find the closest points in the identifiers list for both of the extremes of the entering subpath
       # print "Type of r : " + str(type(r))
        self.add(p,r,Path,i,Id) 

               
    def add(self,p,r,Path,index,ID):
        """
        add the subpath stored in Path[index] to self.Path according to the values present in the tuples p and r , separate 3 cases:\n
        if both are not boolean, which means that both extremes of the subpath have a subpath in self.Path which is closer than the threshold,
        fuse with each of the subpaths one at the time, p[0] represents the subpath to which the left extreme is going to be fused and r[0] the same but for the
        right extreme, r[1] and p[1] refers to the extremes of the closest subpaths to which the extremes of the entering subpath are closest and with it will
        be fused. In the case that r and p ar the same then we just add the one which is closest. Different ways of adding depends in which of the extremes(right
        or left) is going to be added. Similar explanations for the other two cases.
        """
        if type(p)!= bool and type(r)!=bool:

            type2 = setType(p[1])
            if p[0]==r[0] and p[1]==r[1]:
                d1 = ID[0].logdist(self.Identifiers[p[0]][p[1]])
                d2 = ID[1].logdist(self.Identifiers[r[0]][r[1]])
                if d1 >=d2:
                    self.updatePath(Path.subpaths[index],p[0],p[1],ID,'B')
                else:
                    self.updatePath(Path.subpaths[index],r[0],r[1],ID,'F')
            else:
                self.updatePath(Path.subpaths[index],p[0],p[1],ID,'B')  
                self.updatePath(self.subpaths[p[0]],r[0],r[1],self.Identifiers[p[0]],type2)
                
            if p[0] != r[0]:
                self.deleteSubpath(p[0])
            
            
        elif type(p)!= bool or type(r)!=bool:
            if type(p)!= bool:
                self.updatePath(Path.subpaths[index],p[0],p[1],ID,'B')
                
            else:
                self.updatePath(Path.subpaths[index],r[0],r[1],ID,'F')
        else:
            
            self.addnewSubPath(Path,index)
             
                   
    def addnewSubPath(self,Path,index):
        """
        * append the subpath stored in Path[index] to the subpaths list of self.Path
        * append the identifiers(extremes) of the subpath to the identifiers list of self.Path
        """
        self.subpaths.append(Path.subpaths[index])
        
        self.Identifiers.append(Path.Identifiers[index])
                
                        
    def updatePath(self,Update,p,P,ID,Type):
        """
        Update self.Path by adding the subpath Update to the subpath of index p , distinguish two cases
        if the extreme that is going to be fused is the 'rightmost' it is of type 'F'(forward) else is of type
        'B'(backward), it also distinguish cases if the extreme of the subpath to wich Update is going to be fused is
        the 'rightmost' or 'leftmost' , it reverse Update if both extremes are the rightmost(leftmost) 
        """
        if Type == "F":
            if  P == 1:
                self.subpaths[p] += Reverse(Update)
                self.Identifiers[p][1] = ID[0]
            else:
                if self.Identifiers[p] == ID:
                    self.subpaths[p] += [ID[0]]
                else:
                    self.subpaths[p] = Update + self.subpaths[p]
                self.Identifiers[p][0] = ID[0]
            
        else:
            if P == 1:
                if self.Identifiers[p] == ID:
                    self.subpaths[p] += [ID[0]]
                else:
                    self.subpaths[p]+=Update
                self.Identifiers[p][1] =ID[1]
            else:                    
                self.subpaths[p] = Reverse(Update) + self.subpaths[p]
                self.Identifiers[p][0] = ID[1]
        
        
    def getNumSubpaths(self):
        """
        return the number of subpaths
        """
        return len(self.subpaths)
    
    def plot(self,plothandler,linecolor,line,marker):
        SubPaths = [self.getSubPathXY(i) for i in range(len(self.subpaths))]
        for SubP in SubPaths:
            plothandler.plot(np.log10(SubP[0]),SubP[1],linestyle = line ,color = linecolor,marker=marker,linewidth = 1.)
        
#Auxiliar Functions

    
def setType(i):
    """
    Set the type for the fusing method , updatePath
    """
    if i == 0 :
        return 'B'
    else:
        return 'F'

def Reverse(List):
    """ Returns the List with the order of elements reversed"""
    A = [List[i] for i in range(len(List)-1,-1,-1) ]
    return A
def Match(ID,List,Dist):
    """ 
    Find the closest points in List for both elements of ID, if none of the points are at a distance below than Dist, returns a boolean False,
    more precisely if there exist points below the threshold it returns a 2 element tuple whose first element is the index n of List in which the element is located,and the second element is the index within that element of the point for which it has the closest distance.\n
    The function behaves a little different if ID consists of just one point, it search two times for two possible subsets to which join.
    The function assumes that each point in ID(which represents the extremes of a subpath) must join to just one subpath in the case of two different points in ID, and at most to two different subpaths in the other case.
    This assumption is influenced by the fact that we are plotting sign-change boundaries for continuous functions, and intersections of such sets.
    """

    a = ID[0]
    b = ID[1]
    if a!= b:
        return psearch(a,List,Dist),psearch(b,List,Dist)
    else:
        n = psearch(a,List,Dist)
        nList = List[:]
        if type(n)!= bool:
            nList[n[0]][n[1]] = Point(1e-20,-20)
        else:
            return False, False
        n2= psearch(b,nList,Dist)
        return n,n2
        
        

        
def getMinimums(a,List):
    """
    given a Point 'a' and a List of Points returns the minimum distance and the index of List in which it was found.
    """
    index = 0 
    mini = a.logdist(List[0])
    for i in range(1,len(List)):
        d = a.logdist(List[i])
        if d < mini:
            mini = d
            index = i
    return mini,index

def psearch(a,List,Dist):
    """
    from a point a and a List of Two-element list of points, first it calculates the minimum distance to each of the sublists and the element of each sublist
    in which the minimum distance is found and stores both in two new lists dist_list and minIndex_List then finds the min of dist_list, if it is smaller than
    the threshold distance Dist it calculates the index of the min "n" and returns a two element tuple consiting of "n" and the index of the element in the "n" 
    sublist of List for which the minimum was found
    """
    dist_list = []
    minIndex_list = []
    for i in range(len(List)):
        minD,minI = getMinimums(a,List[i])
        dist_list.append(minD)
        minIndex_list.append(minI)
    
    if min(dist_list) < Dist:
        n = min_index(dist_list)
        return (n , minIndex_list[n])
    else:
        return False
    
    
def min_index(dist_list):
    """
    Given a list of numerical entries(distances) , returns the index of the minimum element
    """
    index = 0
    min = dist_list[0]
    for i in range(1,len(dist_list)):
        if ( dist_list[i]) < min:
            min = dist_list[i]
            index = i
    return index
