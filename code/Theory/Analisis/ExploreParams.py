import os
from InvasionAnalysis import *

#######################Explore Parameter Space #####################################

def CreateDirCoder(ParamsToExplore,dimDict):
    """
    ParamsToExplore is a dictionary with parameter names as keys and a list of values as values
    dimDict 
    
    """
    Dict={}
    for K in ParamsToExplore.keys():
        if K == 'fmC' or K == 'fmPC' or K =='fmPR':
            Dict[K]={ParamsToExplore[K][i] : ParamsToExplore[K][i] for i in range(len(ParamsToExplore[K]))}
            
        elif K == 'D_R' or K == 'D_C':

            Dict[K] ={ParamsToExplore[K][i]:dimDict[ParamsToExplore[K][i]] for i in range(len(ParamsToExplore[K]))}
        else:
            Dict[K]={ ParamsToExplore[K][i]: "%.e"%ParamsToExplore[K][i] for i in range(len(ParamsToExplore[K]))}
        
    return Dict

def MakeTotalParamsCombination(ParamsToExplore,AuxiliarParams,SpecialParams):
    """
    
    * ParamsToExplore is a dictionary with parameter names as keys and a list of values as Values
    * AuxiliarParams is a dictionary that stores the information about correlated parameters, that is we set K to a value K0 we must also change the values of the params specified in this dict to the respective values.
    * SpecialParams refer to a superClass of Params which are defined as combinations of the elementary params , e.g Comb1:= {'fmC':Grazing,'fmPC':Active,'fmPR':Grazing}
    
    returns a list of Dictionaries containing all possible combinantions for the given params
    """
    Combinations = [{}]
    
    for param in ParamsToExplore.keys():
        Combinations = AddToCombinations(param,ParamsToExplore,AuxiliarParams,SpecialParams,Combinations)
    
    return Combinations

def AddToCombinations(param,ParamsToExplore,AuxiliarParams,SpecialParams,Combinations):
    """
    From a base set List of Combinations, update it using the values stored in ParamsToExplore[param], AuxiliarParams[param]
    """

    ParamVals = ParamsToExplore[param]
    newCombinations=[]
    for combination in Combinations:
        ValComb = []
        for val in ParamVals:
            Comb = combination.copy()
            Comb = AddValtoComb(Comb,val,param,AuxiliarParams,SpecialParams)
            ValComb+=Comb
        newCombinations+=ValComb
    return newCombinations

def AddValtoComb(Comb,val,param,AuxiliarParams,SpecialParams):
    """
    Find all the associate parameter combinations given a particular value val of the parameter 'param'.
    This is donde in a recursive way, the depth of the recursion is finite since the number of keys for AuxiliarParams is finite.
    If param is a SpecialParam, it adds the low level representation of it to Comb.
    Finally it updates each of the combinations found with Comb and returns it.
    """
    Comb[param] = val
    if param in AuxiliarParams.keys():
        combinations = MakeTotalParamsCombination(AuxiliarParams[param][val],AuxiliarParams,SpecialParams)
    else:
        combinations =[]

    if param in SpecialParams:
        Pars = SpecialParams[param][val]
        for par in Pars :
            Comb[par] = Pars[par]
    return UpdateComb(Comb,combinations)


def UpdateComb(Comb,combinations):
    """
    Update each dictiondary in *combinations* with the dictionary *Comb*.
    """


    if len(combinations)>=1:
        for combination in combinations:
            combination.update(Comb)
        UpdatedCombinations = combinations
    else:
        UpdatedCombinations = [Comb]
        
    return UpdatedCombinations
        
                
    
def ExploreParamSpace(InitDict,ParamsToExplore,TotalParams,xlims,mode,HeaderInv,AuxiliarParams,SpecialParams,dimDict,initDirection,
                      ksim=True,massVals=[0]):
    """
    Explore the parameter space, by first constructing the total params combinations to explore , creating a Dict that contains a name for each combination
    and then Evaluating and Saving each of the parameter combinations
    """
    
    ParamCombinations= MakeTotalParamsCombination(ParamsToExplore,AuxiliarParams,SpecialParams)
    ParamsDirCoder = CreateDirCoder(TotalParams,dimDict)
    CreateTxtCoder(TotalParams,initDirection+"ParamsExplored.txt")
    for mass in massVals:
        for combination in ParamCombinations:
            InitDict.update(combination)
            EvaluateParams(InitDict,mode,xlims,HeaderInv,ParamsDirCoder,initDirection,ksim,mass)


def EvaluateParams(paramdict,mode,xlims,HeaderInv,ParamsDirCoder,Direction,ksim,mass):
    r"""
    For a given dict of parameters:
    * calculate the invasibility boundaries, the Equilibrium, eigenvalues ,MTP, Zones and width of the zones for (x,y) values between :math:`xRange \times [-10,5`
    * Write the results to csv files
    """
    paramdict['massP'] = mass
    #Construct Object
    WD = BSR(paramdict,mode,xlims,ksim)
    WD.setfDict()
    #Set Focus 
    if ( ksim == False and paramdict['fmC'] == 'Grazing' ):
        WD.getandSetxFocus(300,5e-4)
    
    Inv = InvBoundaries(WD)   
    if ( ksim == False and paramdict['fmC'] == 'Grazing' ):
        Inv.setUpGuess(15)
      
    #create directions
    Types = ['Inv','Widths','Zones']
    
    DirInv,DirWidths,DirZones = [Direction+ConstructDir(paramdict,mode,ParamsDirCoder,Type=T) for T in Types]
    
    #Invasibility
    Inv.setAndWriteInvBoundaries(HeaderInv,DirInv) 
    Inv.setPositiveBoundaries()
    Inv.setAndWriteWidthsZones(DirWidths,DirZones)

     
    
def ConstructDir(ParamDict,mode,ParamsDirCoder,Type):
    """
    Construct a Dictionary of file directions using the values of the parameters in ParamDict and the map in ParamsDirCoder
    """
    Dir = ""
    sufix=".csv"
    massCode = "massP"+"%.e"%(ParamDict['massP'])
    for param in ParamsDirCoder.keys():
        u = param+str(ParamsDirCoder[param][ParamDict[param]])
        Dir+=u

    return Type+mode+Dir+massCode+sufix 

    


def CreateTxtCoder(ParamsDirCoder,direction):
    """
    Create a text file in the specified direction, with the details related to all the parameters used in the computation.
    """
    F = open(direction,'w')
    for Key in ParamsDirCoder.keys():
        Line = str(Key) + ' : '
        for val in ParamsDirCoder[Key]:
            Line += str(val) + '\t'
        Line+='\n'
        F.write(Line)
    F.close()


