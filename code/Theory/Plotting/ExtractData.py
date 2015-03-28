import os
import re

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


def getDataSets(Directory,params,Types,Dim,masses,TypeDict):
    """
    From a Directory of files,even other folders, for each of the elements in Dim(which will usually denote dimensions, but it could also be any other focal parameter,extract using the ExtractDataSets method : for each of the mass values in masses , all the datasets names which match the parameter values given in params.
    """
    DataSets=[]
    for D in Dim:
        param=params.update(D)
        Data=[]
        for mass in masses:
            param.update(mass)
            RawData= ExtractDataSets(Directory,param,Types)
            Data.append(FormatData(RawData,TypeDict))
        DataSets.append(Data)
    return DataSets


def GetDataSets(Directory,params,Types,p1,p2,p3,TypeDict):
    """
    p_i is a list of dicts, params is the baseline parameter values used in the search
    """
    DataSets=[]
    for p in p3:
        params.update(p)
        Data = getDataSets(Directory,params,Types,p2,p1,TypeDict)
        DataSets.append(Data)
    return DataSets


def FormatData(RawData,TypeDict):
    NewData ={}
    for Type in RawData:
        D =TypeDict[Type](RawData[Type][0])
        NewData[Type]=D
    return NewData


        

