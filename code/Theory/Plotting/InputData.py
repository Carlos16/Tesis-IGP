import csv
 
class Data(object):
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
        paramsEspecifications = data_handler[-1]
        distribution = data_handler[-2]
        output = InvData(dataset,InvScenarios,paramsEspecifications,distribution)
        
        return output
         
       
        
class InvData(Data):   
    def __init__(self,data,InvScenarios,paramsEspecifications,distribution):
        Data.__init__(self,data)
        self.Scenarios = InvScenarios
        self.pEspecifications = paramsEspecifications    
        self.distribution = distribution
        self.formated_data = ''
        self.OrderedData = {}
        self.Paths= {}
    def formatData(self):
        """
        Transform the original data to a format amenable for future working, first transform it to a 2 Dimensional array
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
    def PlotBoundaries(self,plothandler):
        self.OrderData()
        self.Plot()
        
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
        YPath.FormatPoints()
        return YPath
        
            
        
        
    
    def Order(self,scenario):
        Set= self.formated_data[scenario]
        new_x = []
        new_y = []
        if len(Set['x'])>1:
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
        if d>0.55:
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
    
