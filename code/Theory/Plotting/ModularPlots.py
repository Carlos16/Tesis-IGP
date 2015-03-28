from PointsPath import *
from ReadData import *
from ExtractData import *

from numpy import ndarray

def LoadInv(dataset):
    P = InputInvData(dataset)
    return P.Read(',')
def LoadMTP(dataset):
    P = InputMatrixData(dataset)
    return P.Read(',',MTPData)

def LoadStab(dataset):
    P = InputMatrixData(dataset)
    return P.Read(',',EigenValData)

def MakePlots(ksim,Types,params,Directory, Dim = ['2D','3D'],masses=['mass1e-10','mass1e-05','mass1e+00','mass1e+05'],
              TypeDict={'Inv':LoadInv,'MTP':LoadMTP,'Stab':LoadStab},figsize=(10,8),
              paramDict={'Inv':[ColorCoder,LineCoder],'MTP':[yRange,xSep,ySep,mapcolor],'Stab':[yRange,xSep,ySep,mapcolor]}):
    
    fig,axes = makeFigAndAxes(ksim,figsize)
    paramDict = formatParamD(paramDict,ksim)
    DataSets = getDataSets(Directory,params,Types,Dim,masses,TypeDict)
    PlotDataSets(axes,DataSets,ksim,paramDict)
    

def formatParamD(paramDict,ksim):
    newParamDict = paramDict.copy()
    for Type in paramDict:
        if not(Type == 'Inv'):
            if ksim:
                newParamDict[Type] = [[-14,6]]+paramDict[Type]
            else:
                newParamDict[Type] = [[-10,5]]+paramDict[Type]
    return newParamDict
                    
def makeFigAndAxes(ksim,Figsize):
    if ksim :
        return plt.subplots(nrows=len(masses),ncols=len(Dim),sharex=True,sharey=True,figsize=Figsize)
    else:
        return plt.subplots(nrows=len(Dim),ncols=len(masses),sharex=True,sharey=True,figsize=Figsize)
    
    
def PlotDataSets(axes,DataSets,ksim,paramDict):
    for i in range(len(axes)):
        PlotD(axes[i],DataSets[i],ksim,paramDict)

        
def PlotD(ax,Data,ksim,paramDict):
    if type(ax) !=  ndarray:
        Plot(ax,Data[0],paramDict)
    else:
        for j in range(len(ax)):
            Plot(ax[j],Data[j],paramDict)
            
def Plot(ax,dataset,paramDict):
    for Type in dataset:
        if Type == 'Stab':
            dataset['Stab'].plotStabBound(ax,*paramDict['Stab'])
        else:
            dataset[Type].plot(ax,*paramDict[Type])
        
               

