
class MyFloat(float):
    """differs from float in just the printing method"""
    def __str__(self):
        if self == 0.0:
            return "0"
        else:
            return "%.20f"%self

class MyTuple(tuple):
    """differs from tuple in just the printing method"""
    def __str__(self):
        string = str(self[0])
        for i in range(1,len(self)):
            string = string + ':'+str(self[i])
          
        return '('+string+')'
class MyInnerTuple(tuple):
    def __str__(self):
        string = str(self[0])
        for i in range(1,len(self)):
            string = string + '|'+str(self[i])
          
        return '('+string+')'
