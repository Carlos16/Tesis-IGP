import csv
import PointsPath
class Data(object):
    def __init__(self,data,paramsEspecifications,xFocus,xFSep):
        self.data = data
        self.paramsEspecificiations = paramsEspecifications
        self.xFocus = TransformToFloat(xFocus,':')
        self.xFSep = float(xFSep)
     
          
        
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
        
        
    def plot(self,plothandler,xlims,ylims,xsep,ysep,mapcolor):
        self.formatData()
        #create matrix
        X = setXRange(xlims,self.xFocus,self.xFSep,xsep)
        Y = np.arange(ylims[0],ylims[1],ysep)
        X,Y = np.meshgrid(X,Y)
        #create levels and norm
        lev,norm1 = self.setLevelsandNorm()
        #plot heatmap
        plothandler.contourf(X,Y,self.getPlotObject(),cmap = mapcolor,levels=lev,norm=norm1)
        plothandler.contourf(X,Y,self.getPlotObject(),cmap = mapcolor,levels=lev,norm=norm1)

class InputInvData(Data):
    def __init__(self,direction):
        self.direction = direction
    def Read(self,delimiter):
        with open(self.direction,'rb') as csvfile:
            reader = csv.reader(csvfile,delimiter=',',quotechar='|')
            data_handler = []
            for row in reader:
                data_handler.append(row)
        dataset = np.array(data_handler[1:-2],dtype=object)
        InvScenarios = data_handler[0]
        paramsEspecifications = data_handler[-2]
        distribution = data_handler[-3]
        output = InvData(dataset,InvScenarios,paramsEspecifications,distribution)
        return output

class InputMatrixData(Data):
    def __init__(self,direction):
        self.direction = direction
        
    def Read(self,delimiter,dataType):
        with open(self.direction,'rb') as csvfile:
            reader = csv.reader(csvfile,delimiter=',',quotechar='|')
            data_handler = []
            for row in reader:
                data_handler.append(row)
        dataset = np.array(data_handler[0:-2],dtype=object)
        dataset.reshape((len(dataset),len(dataset[0])))
        paramsEspecifications = data_handler[-2]
        xFocus,xFSep = data_handler[-1]
        
        output = dataType(dataset,paramsEspecifications,xFocus,xFSep)
        
        return output

    
    
def TransformToFloat(X,delimiter):
    return [ float(item) for item in X[1:-1].split(delimiter)]

    
class EqData(Data):
    def __init__(self,data,paramsEspecifications,xFocus,xSep):
        Data.__init__(self,data,paramsEspecifications,xFocus,xSep)
        self.formated_data = ''        
        
    def extract_Data(self):
        self.TransformtoFloats()
        m,n = self.data.shape
        R = np.zeros((m,n))
        C = np.zeros((m,n))
        P = np.zeros((m,n))
        for i in range(m):
            for j in range(n):
                X = self.data[i,j]
                R[i,j]=X[0]
                C[i,j]=X[1]
                P[i,j]=X[2]
        return R,C,P

    
class EigenValData(Data):
    def __init__(self,data,paramsEspecifications,xFocus,xSep):
        Data.__init__(self,data,paramsEspecifications,xFocus,xSep)
        self.DominantEigenVal = 0.
        self.numEigenVal = 0
        self.StabBound =0
        
    def TransformtoFloats(self):
        m,n = self.data.shape
        for i in range(m):
            for j in range(n):        
                self.data[i,j] = [ TransformToFloat(item,'|') for item in self.data[i,j][1:-1].split(':')]
                
    def formatData(self):
        self.TransformtoFloats()
        self.setDominants()
        self.setStabBound()
        
    def setStabBound(self):
        self.StabBound = getStabBoundary(self.DominantEigenVal)
        
    def plotStabBound(self,plothandler,xlims,ylims,xsep,ysep,mapcolor):
        self.formatData()
        #create matrix
        X = setXRange(xlims,self.xFocus,self.xFSep,xsep)
        Y = np.arange(ylims[0],ylims[1],ysep)
        X,Y = np.meshgrid(X,Y)
        #create levels and norm
        lev= [-0.5,0.5,1.5]
        #plot heatmap
        plothandler.contourf(X,Y,self.StabBound,cmap = mapcolor,levels=lev)

    def getPlotObject(self):
        return self.DominantEigenVal
    def numEigs(self):
        m,n = self.data.shape
        Num = np.zeros((m,n))
        for i in range(m):
            for j in range(n):
                Num[i,j]=len(self.data[i,j])
                
        self.numEigenVal = Num

    def setDominants(self):
        m,n = self.data.shape
        Dom = np.zeros((m,n))
        for i in range(1,m):
            for j in range(n):
                k = getMaxRealIndex(self.data[i,j])
                
                Dom[i,j]=self.data[i,j][k][0]
        for j in range(n):
            Dom[0,j] = Dom[1,j]
                
        self.DominantEigenVal = Dom
        
    
        
    def setLevelsandNorm(self):
        m,M = getMinMax(self.DominantEigenVal)
        lev = setLevel(m,M,100)
        norm = mpl.colors.BoundaryNorm(lev,256)
        return lev,norm
    
    
    
    
def setLevel(m,M,npoints,zero_=1e-18):
    if M >0:
        if m <0 :
            if abs(m)>zero_:
                lev1 = -(10**np.linspace(np.log10(zero_),np.log10(abs(m)),npoints/2))[::-1]
                lev2 =  10**np.linspace(np.log10(zero_),np.log10(M),npoints/2)
                return np.concatenate([lev1,lev2])
            else:
                return 10**np.linspace(np.log10(zero_),np.log10(M),npoints)
        else:
            return 10**np.linspace(np.log10(m),np.log10(M),npoints)
                
        
    else:
        return -(10**np.linspace(np.log10(abs(M)),np.log10(abs(m)),npoints))[::-1]
        
        
class MTPData(Data):
    def __init__(self,data,paramsEspecifications,xFocus,xSep):
        Data.__init__(self,data,paramsEspecifications,xFocus,xSep)
        
    def TransformtoFloats(self):
        m,n = self.data.shape
        for i in range(m):
            for j in range(n):
                
                self.data[i,j] = float(self.data[i,j])
                
    def formatData(self):
        self.TransformtoFloats()
        
    def getPlotObject(self):
        return self.data
       
    def setLevelsandNorm(self):
        L = np.arange(2.1,3.05,0.05)
        lev = [0.,1.,2.] +list(L)
        norm1 = mpl.colors.BoundaryNorm(lev,256)
        return lev,norm1
        
        
                      
def setXRange(xlims,xfocus,xfocussep,xsep):
    A = np.arange(xlims[0],xfocus[0],xsep)
    B = np.arange(xfocus[0],xfocus[1],xfocussep)
    C = np.arange(xfocus[1],xlims[1],xsep)
    D = np.concatenate([A,B,C])
    return D

def getMaxRealIndex(List):
    Max = List[0][0]
    index = 0 
    for i in range(1,len(List)):
        if(List[i][0]>Max):
            index = i
            Max = List[i][0]
    return index             
        
        
    
    
def isitCoexistence(X,index,index2):
    try:
        return isinPositiveOrthant(X,index)
    except:
        try:
            return isinPositiveOrthant(X,index2)
        except:
            return isinPositiveOrthant(X,0)
        
    
def isitPositive(X):
    pos = 0
    for i in range(len(X)):
        p = isinPositiveOrthant(X,i)
        if p>0:
            return p
    return pos

def isinPositiveOrthant(X,index):
    P = X[index][2]
    if P>0:
        return 1
    else:
        return 0
            
        
         
class InvData(Data):   
    def __init__(self,data,InvScenarios,paramsEspecifications,distribution):
        Data.__init__(self,data,paramsEspecifications,'(0:0)',0.)
        self.Scenarios = InvScenarios
        self.distribution = distribution
        self.formated_data = ''
        self.OrderedData = {}
        self.Paths= {}
        
    def formatData(self):
        """" Transform the original data to a format amenable for future working, first transform it to a 2 Dimensional array
        whose number of columns if the number of distinct scenarios encountered in the invasibility analysis(which is four in our case)
        and the number of rows is the maximum number of points present in the invasibility set of any of the scenarios. It also convert
        the string elements to floats, and returns a dictionary of dictionaries for each of the scenarios with x and y coordinates as
        the two keys of each of the dictionaries and whose elements are a list of lists for each of the distinct invasibility sets 
        contained within each scenario. It furthers filters the results for the scenarios P to C-R and C to P-R neglecting all the
        elements in which the necessary conditions,C to R and P to R respectively, are fullfiled."""
        self.reshape()
        self.TransformtoFloats()
        self.distribution = [ [int(item) for item in self.distribution[i][1:-1].split(':')] for i in range(len(self.distribution)) ]
        distribution = self.distribution
        formated_data = {}
        ncols = self.ncols()
        for n in range(ncols):
            colDist = distribution[n]
            formated_data[self.Scenarios[n]] = {'x':[],'y':[]}
            for i in range(len(colDist)):
                x_data =[]
                y_data =[]
                for j in range(sum(colDist[0:i]),sum(colDist[0:i+1])):
                    x_data.append(self.data[j,n][0])
                    y_data.append(self.data[j,n][1])
                formated_data[self.Scenarios[n]]['x'].append(x_data)
                formated_data[self.Scenarios[n]]['y'].append(y_data)
                        
        
        self.formated_data = formated_data
    def plot(self,plothandler,colorCoder,lineCoder):
        self.formatData()
        self.OrderData()
        self.InnerPlot(plothandler,colorCoder,lineCoder)
        
        
    def InnerPlot(self,plothandler,colorCoder,lineCoder):
        for key in self.Paths:
            self.Paths[key].plot(plothandler,colorCoder[key],lineCoder[key])
         
             
    def FilterData(self,Target,xFilter,yFilter):
        """Filter the data contained in the target key buy deleting all the number which are below
        the yFilter correspondence of each of the xFilter positions, this implimentation use a
        binary search algorithm for finding the elements of xFilter within each of the sublists storaged
        in Target 
        @param Target  a key of the data dictionary
        @param xFilter a list of x coordinates
        @param yFilter a list of y coordinates
        """
        
        T = self.formated_data[Target]
        Xlists = T['x']
        Ylists = T['y']
        
        for x_i in range(len(xFilter)):
            x = xFilter[x_i]
            
            for index in range(len(Xlists)):
                    
                p = BSearch(Xlists[index],x)
                if type(p)!=bool:
                    if yFilter[x_i] < Ylists[index][p]:
                        break
                    else:
                        Ylists[index].pop(p)
                        Xlists[index].pop(p)
                else:
                    break   
    
    def OrderData(self):
        Data = self.formated_data
        for scenario in Data.keys():
            self.Order(scenario)
        self.ConstructPaths()
        
    def ConstructPaths(self):
        
        for scenario in self.Scenarios:
            Path = self.BuildPath(scenario)
            self.Paths[scenario] = Path
        
    def BuildPath(self,scenario):
        P = self.OrderedData[scenario]
        YPath = Path(P['x'][0],P['y'][0])
        for index in range(len(P['x']) -1 ) :
            YPath.AddPath(Path(P['x'][index+1],P['y'][index+1]))
        #YPath.FormatPoints()
        return YPath
        
          
    def Order(self,scenario):
        Set= self.formated_data[scenario]
        new_x = []
        new_y = []
        if len(Set['x'])>=1:
            for index in range(len(Set['x'])):
                
                new_subx,new_suby= self.Classify(scenario,index)
                new_x.append(new_subx)
                new_y.append(new_suby)
        else:
            new_x.append([Set['x'][0]])
            new_y.append([Set['y'][0]])
        self.OrderedData[scenario] = {'x':new_x,'y':new_y}
        
    
        
    
    def Classify(self,scenario,index):
        
        xset = self.formated_data[scenario]['x'][index]
        yset = self.formated_data[scenario]['y'][index]
        
        #searchRange = self.getIndexes(scenario,index)
        
        newX = []
        newY = []
        roamingX=[xset[0]]
        roamingY=[yset[0]]
        for i in range(0,len(xset)-1):
            if self.isDist(xset,yset,i):
            #if self.isDist(xset,yset,searchRange,i):
                roamingX.append(xset[i+1])
                roamingY.append(yset[i+1])
            else:
                newX.append(roamingX)
                newY.append(roamingY)
                roamingX=[xset[i+1]]
                roamingY=[yset[i+1]]
            
        newX.append(roamingX)
        newY.append(roamingY)
        
        return newX,newY
    
    def isDist(self,xset,yset,i):
        P0 =  Point(np.log10(xset[i]),yset[i])
        P1 = Point(np.log10(xset[i+1]),yset[i+1])
        d = P0.dist(P1)
        if d>0.2:
            return False
        else:
            return True
        
        
    
    def getIndexes(self,scenario,index):
        n = len(self.formated_data[scenario]['x'])
        t = range(n)
        t.pop(index)

        
        searchRange = {}
        for i in  t:
            searchRange[i]={'x':self.formated_data[scenario]['x'][i],'y':self.formated_data[scenario]['y'][i]}
        return searchRange
   
