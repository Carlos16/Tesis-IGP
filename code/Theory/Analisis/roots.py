from scipy.optimize import brentq
import numpy as np

###Code for computing roots of the function using the
### brentq method
def Get_roots(f,arg,array,j=0,debug=False,method=brentq,xtol=1e-30):
    """
    Search for the all the zeros of a continuos and smooth function f, using the brentq optimization algorithm
    * Assumes that there is at most one zero in any of the subintervals (x[i],x[i+1]), where d(x[i+1],x[i]) <=0.03 
    
    The way it performs the search is the following:
    1) Computes F_array using the function f, array and any other arguments that the functions needs(This means that f could be a multivariate function,
    but we are just fixing a subset of a line in which to search), this arguments are formated by the createArray function and we are making use of the convenient
    way python functions handle numpy arrays.
    2) set a inital index j to 0
    3) Starts at F_array[j] computes the sign of it, proceed through the elements of F_array until it founds i such that sign(F_array[j]) != sign(F_array[i]) 
    4) apply the brentq method to the interval array[j],array[i] and give all the necessary arguments.
    5) changes j to i
    6) Repeats steps 3-4-5 until there are no more elements in F_array
    """
    n=len(array)
    if debug:
        print len(array)

    F_array=f(array,*createArray(arg,array))

    if debug:
        print len(F_array)
    roots=[]
    while j< n -1:
        sign_ = np.sign(F_array[j])
        for i in range(j,n):
            if np.sign(F_array[i])!= sign_:
                roots.append(np.log10(method(f,array[j],array[i],args=tuple(arg),xtol=xtol)))
                j = i 
                break
            if i == n-1:
                j = i
    return roots
def createArray(args,baseArray):
    """
    Create a list of np.arrays whose values are all the same, using the values stored in args
    """
    n = len(baseArray)
    return [np.array([args[i]]*n) for i in range(len(args))]
