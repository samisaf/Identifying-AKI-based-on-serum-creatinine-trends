from pandas import Series, DataFrame
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

ptDict = {} # Dictionary containing all lab values
Patients = {} # Dictonary, keys are MRN, values are Patient objects

class Patient(object):
    def __init__(self, mrn: int, esrd: bool = False,  crs: dict = dict()):
        self.mrn = mrn
        self.esrd = esrd
        self.crs = DataFrame(Series(crs), columns=['value'])
        self.crs.index = pd.to_datetime(self.crs.index)
        self.baseCr = self.__calcBaseCr__()
        self.crs['order'] = range(len(crs))
        self.crs['slope'] = self.__calcSlopes__()
        self.crs['peak'] = self.__calcPeaks__()
        self.crs['aki'] = self.__calcAKI__()
        self.aki = self.crs.value[self.crs.aki]
 
    def __str__(self): 
        temp = "<MRN {}, baseline cr {}, crs {}, AKI {}>"\
            .format(self.mrn, self.baseCr, len(self.crs), sum(self.crs.aki))
        return temp

    def __repr__(self): return self.__str__()
    
    def __calcBaseCr__(self):
        return np.percentile(self.crs.value, 50)
    
    def __calcSlopes__(self):
        slopes = []
        for i in range(len(self.crs)-1):
            x1, x2 = self.crs.order[i], self.crs.order[i+1]
            y1, y2 = self.crs.value[i], self.crs.value[i+1]
            slopes.append((y2-y1)/(x2-x1))
        # assign slope after last point to the preceding slope
        if len(slopes) == 0: slopes.append(0) 
        else: slopes.append(slopes[len(slopes) - 1])
        return slopes
        
    def __calcPeaks__(self):
        # first point is a peak if following slope is negative
        peaks = [self.crs.slope[0] <= 0] 
        for i in range(1, len(self.crs)):
            if self.crs.slope[i-1] > 0 and self.crs.slope[i] <= 0: peaks.append(True)
            else: peaks.append(False)
        return peaks
    
    def __calcAKI__(self):
        return (self.crs.value > self.baseCr * 1.5) & self.crs.peak
                               
    def plot(self):
        plt.hlines(self.baseCr, self.crs.index.min(), self.crs.index.max(), linestyles='dotted')
        plt.plot(self.crs.index, self.crs.value)        
        plt.plot(self.crs.index, self.crs.value, 'g.')
        plt.plot(self.crs.index[self.crs.aki], self.crs.value[self.crs.aki], 'ro')
        plt.title(str(self))
        plt.xlabel("Date")
        plt.ylabel("Creatinine")
        plt.show()
     
def createCrDict(df: DataFrame):
    global ptDict 
    uniqueMRN = np.unique(df.MRN)
    for i in uniqueMRN: 
        ptDict[i] = dict()
    col = 1
    while col < len(df.columns):
        creatinines = df.ix[:, col]
        dates = df.ix[:, col + 1]
        for i in range(0, len(df)):
            mrn = df.MRN[i]
            c = creatinines[i]
            d = dates[i]
            if not(np.isnan(c)): ptDict[mrn][d] = c
        col = col + 2
    return ptDict

def createPts(patients: dict):
    global Patients
    for key in patients:
        temp = Patient(mrn=key, crs=patients[key])
        Patients[key]=(temp)

def plotPatients():
    global Patients
    for key in Patients:
        p = Patients[key]
        plt.figure()
        p.plot()

def getNumAKI():
    global Patients
    mrn = list()
    aki = list()
    for key in Patients:
        mrn.append(Patients[key].mrn)
        aki.append(Patients[key].aki.size)
    di = {'MRN': mrn, 'NumAKI': aki}
    df = DataFrame(di)
    df['AKI'] = df.NumAKI > 0
    return DataFrame(df)

def start(path:str):
    global Patients  
    df = pd.read_csv(path)
    createCrDict(df)
    createPts(ptDict)
    getNumAKI()

if __name__ == "__main__":
    if len(sys.argv) > 1: start(sys.argv[1])
    else: 
        path = os.getcwd() + "\\" + "Labs.csv"
        start(path)
    input("Done, press enter to exit...")
