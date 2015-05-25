from ExploreParams import *

"""
This module explore the parameter space

"""

#Base Dictionary#

params={'w':0.75,'pd':0.21,'pv':0.26,'Er':0.,'Ek':0.,'ER':0.,
        'EC':0.,'EP':0.,'Eq1':0.,'Eq2':0.,'TR':300.,'TC':300,
        'TP':300.,'D_R':2.,'D_C':2.,'fmC':'Grazing','thermyR':'Ecto',
        'thermyC':'Ecto','thermyP':'Ecto','fmPC':'Active','fmPR':'Grazing','k0':0.01,'r0':1.71e-6,'a012':10**(-3.08),
        'a03':10**(-3.08),'d0':3.,'q10':4.15e-8,'q20':4.15e-8,'v0R':1.,'v0C':1,'v0P':1,'k':b_k,'e1':0.3,'e2':0.1,'e3':0.5,
        'hC0':1.,'hP0':1.,'formPC':2,'formPR':1,'formC':1,'a':1.,'b':2.,'c':0}


#Exploration#

ParamsToExplore= { 'D':['2D','3D'],'fm':['comb1','comb2','comb3']}
dimDict = {2:'2D',3:'3D'}
TotalParamsExplored = { 'D_C':[2.,3.],'D_R':[2.,3.],'k0':[3.,30.,300.,0.01,0.1,1.],'b':[0.02,0.2,2.],'fmPC':['Active','Sit'],'fmC':['Active','Grazing'],'fmPR':['Active','Grazing','Sit'],'e1':[0.3,0.6],'e2':[0.1,0.4],
                      'e3':[0.5]}
ParamsDirCoder = CreateDirCoder(TotalParamsExplored,dimDict)
AuxiliarParams = {'D':{'2D':{'k0':[0.01,0.1,1.],'a012':[10**(-3.08)],'a03':[10**(-3.08)],'pd':[0.21]},
                       '3D':{'k0':[3.,30.,300.],'a012':[10**(-1.77)],'a03':[10**(-1.77)],'pd':[0.2]}},
                  'fm':{'comb1':{'fbs':['comb1'],'e':['comb1']},'comb2':{'fbs':['comb2'],'e':['comb23']},'comb3':{'fbs':['comb3'],'e':['comb23']}},
                  'fbs':{'comb1':{'b':[0.02,0.2,2.]},'comb2':{'b':[0.02,0.2,2.]},'comb3':{'b':[0.02,0.2,2.]}}}
SpecialParams= {'D':{'2D':{'D_C':2.,'D_R':2.},'3D':{'D_C':3,'D_R':3.}},
                'fm':{'comb1':{'fmPC':'Active','fmC':'Grazing','fmPR':'Grazing'},
                      'comb2':{'fmPC':'Active','fmC':'Active','fmPR':'Active'},
                      'comb3':{'fmPC':'Sit','fmC':'Active','fmPR':'Sit'}},
                'fbs':{'comb1':{'formPC':2,'formPR':2,'formC':2},
                       'comb2':{'formPC':2,'formPR':2,'formC':2},
                       'comb3':{'formPC':2,'formPR':2,'formC':2}
                      },
               'e':{'comb1':{'e1':0.3,'e2':0.1,'e3':0.5},
                    'comb23':{'e1':0.6,'e2':0.4,'e3':0.5}}}
Comb = MakeTotalParamsCombination(ParamsToExplore,AuxiliarParams,SpecialParams)
initDirection="c:/Users/Carlos/Documents/Tesis-IGP/Data/InvariantSR/"
 
xlims = [-13,7]
mode = 'LV'
Header= ['Inv C2','Inv P3','Inv P4','Inv C5']
ksim = True
massVals = [1e5]
if __name__== '__main__':
    ExploreParamSpace(params,ParamsToExplore,TotalParamsExplored,xlims,mode,Header,AuxiliarParams,SpecialParams,dimDict,initDirection,ksim,massVals)





