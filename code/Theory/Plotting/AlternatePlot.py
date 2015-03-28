
"""
File : AlternatePLot.py
-----------------------
This file construct a series of plots arranged in a matrix format, whose size is specified by the length of p2(rows) and p3(columns) respectively, p1 determines the different curves in each plot. It first search and arrange the data sets matching the parameter values specified by p1,p2,p3 and the basaline params in a list of lists of lists of dictionarys.Then it creates each plot in a iterative manner, this code assumes that the Matrix of plots has at least dimension 4.
"""




def PlotInt(Directory,params,Types,p1,p2,p3,TypeDict,sizePlot,plotParamsDict,Zones):
    fig,axes = makeFigure(p2,p3)
    DataSets = GetDataSets(Directory,params,Types,p1,p2,p3,TypeDict)
    plotData(Datasets,axes,plotParamsDict,Zones)

def makeFigure(p2,p3,sizePlot):
    fig,axes = plt.subplots(nrows=len(p2),ncols=len(p3),sharex =True,sharey=True,figsize = sizePlot)


def plotData(Datasets,axes,plotParamsDict,Zones):
    for i in len(Datasets):
        for j in Dataset[i]:
            InnerPlot(Dataset[i][j],axes[i,j],plotParamsDict)


def InnerPlot(Datasets,canvas,plotParamsDict):
    for Data in DataSets:
        for Type in Types:
            Data[Type].plot(canvas,*plotParamsDict[Type])
        
        

    
