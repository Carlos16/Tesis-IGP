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

