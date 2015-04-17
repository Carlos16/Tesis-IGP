import os
import re
from ReadData import *

def ExtractDataSets(Dir,params,Types):
    """
    From a Directory, extract all the datasets names that were computed using the parameter values given in params, Types refer to the different file formats that are for each set of parameters.
    """
    ListDirs = os.listdir(Dir)
    ListTypes={T:[] for T in Types}
    for param in params:
         ListDirs= ParseParams(ListDirs,param,params[param])
    for Type in Types:
        ListTypes[Type] = ParseParams(ListDirs,Type,[''])
    return ListTypes
    
    
def ParseParams(ListDirs,param,paramvals):
    """
    From a set of files in a Dir, ListDirs, extract all the files names that were computed using for the parameter param , the values specified in paramvals.
    """
    Items=[]
    for f in ListDirs:
        for val in paramvals:
            match = re.match('.*'+param+str(val),f)
            if match:
                Items.append(f)
    return Items


def getDataSets(Directory,params,Types,p,p2,p1,TypeDict):
    """
    From a Directory of files,even other folders, for each of the elements in Dim(which will usually denote dimensions, but it could also be any other focal parameter)extract using the ExtractDataSets method : for each of the mass values in masses , all the datasets names which match the parameter values given in params.
    """
    DataSets=[]
    p2Name = p2.keys()[0]
    p1Name = p1.keys()[0]
    for par2 in p2[p2Name][p]:
        Data = []
        params[p2Name] = [par2]
        for par1 in p1[p1Name]:
            params[p1Name] = [par1]
            RawData= ExtractDataSets(Directory,params,Types)
            Data.append(FormatData(Directory,RawData,TypeDict))
        DataSets.append(Data)
    return DataSets


def GetDataSets(Directory,params,Types,p1,p2,p3,TypeDict):
    """
    p_i is a list of dicts, params is the baseline parameter values used in the search
    """
    DataSets=[]
    p3Name = p3.keys()[0]
    for p in p3[p3Name]:
        params[p3Name] = [p]
        Data = getDataSets(Directory,params,Types,p,p2,p1,TypeDict)
        DataSets.append(Data)
    return DataSets

def FormatData(Dir,RawData,TypeDict):
    NewData ={}
    for Type in RawData:
        D =TypeDict[Type](Dir+RawData[Type][0])
        NewData[Type]=D
    return NewData


p3 ={'D_C':['2D','3D']}
p2 ={'k0':{'2D':['1e-02','1e-01',r'1e\+00'],'3D':[r'3e\+00',r'3e\+01',r'3e\+02']}}
p1 ={'massP':['1e-10','1e-05',r'1e\+00',r'1e\+05']}

Dir ='C:/Users/Carlos/Documents/Tesis-IGP/Data/VariantSR/'

params = {'fmPR':['Grazing'],'b':[r'2e-02']}
Types = ['Zones']
if __name__== '__main__':
    D = GetDataSets(Dir,params,['Zones'],p1,p2,p3,{})
    print D



