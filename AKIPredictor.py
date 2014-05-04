from pandas import Series, DataFrame
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob

ptDict = {} # Dictionary containing all lab values
Patients = {} # Dictonary, keys are MRN, values are Patient objects

class Patient(object):
    def __init__(self, mrn: int, data: dict = dict()):
        self.mrn = mrn
        self.data = DataFrame(Series(data), columns=['value'])
        self.data.index = pd.to_datetime(self.data.index)
        self.crs = self.data.value
        self.baseCr = self.__calcBaseCr__()
        self.data['order'] = range(len(data))
        self.data['slope'] = self.__calcSlopes__()
        self.data['peak'] = self.__calcPeaks__()
        self.data['aki'] = self.__calcAKI__()
        self.aki = self.data.value[self.data.aki]
 
    def __str__(self): 
        temp = "<MRN {}, baseline cr {}, crs {}, AKI {}>"\
            .format(self.mrn, self.baseCr, len(self.crs), sum(self.data.aki))
        return temp

    def __repr__(self): return self.__str__()
    
    def __calcBaseCr__(self):
        return np.percentile(self.crs, 25)
    
    def __calcSlopes__(self):
        slopes = []
        for i in range(len(self.crs)-1):
            x1, x2 = self.data.order[i], self.data.order[i+1]
            y1, y2 = self.data.value[i], self.data.value[i+1]
            slopes.append((y2-y1)/(x2-x1))
        # assign slope after last point to the preceding slope
        if len(slopes) == 0: slopes.append(0) 
        else: slopes.append(slopes[len(slopes) - 1])
        return slopes
        
    def __calcPeaks__(self):
        # first point is a peak if following slope is negative
        peaks = [self.data.slope[0] <= 0] 
        for i in range(1, len(self.crs)):
            if self.data.slope[i-1] > 0 and self.data.slope[i] <= 0: 
                peaks.append(True)
            else: peaks.append(False)
        return peaks
    
    def __calcAKI__(self):
        return (self.data.value > self.baseCr * 1.5) & self.data.peak
                               
    def plot(self):
        plt.hlines(self.baseCr, self.data.index.min(), self.data.index.max(), linestyles='dotted')
        plt.plot(self.data.index, self.data.value)        
        plt.plot(self.data.index, self.data.value, 'g.')
        plt.plot(self.data.index[self.data.aki], self.data.value[self.data.aki], 'ro')
        plt.title(str(self))
        plt.xlabel("Date")
        plt.ylabel("Creatinine")
        plt.show()
        
def getNumPatients():
    return len(Patients)
    
def getNumCrs():
    return sum([Patients[p].crs.size for p in Patients])
     
def createPtDict(df: DataFrame):
    global ptDict 
    uniqueMRN = np.unique(df.ix[:, 0])
    for i in uniqueMRN: ptDict[i] = dict()
    for i in range(len(df)):
        mrn = df.ix[i, 0]
        creatinine = df.ix[i, 1]
        date = df.ix[i, 2]
        if not(np.isnan(creatinine)): ptDict[mrn][date] = creatinine
    return ptDict

def createPts(patients: dict):
    global Patients
    for key in patients:
        crs = patients[key]
        mrn = int(key)
        if len(crs) > 0: Patients[mrn] = Patient(mrn=mrn, data=crs)

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
    di = {'MRN': mrn, 'NumOfAKI': aki}
    df = DataFrame(di)
    df['AKI'] = df.NumAKI > 0
    return DataFrame(df)

def read(files:[str]):
    global Patients
    for file in files: createPtDict(pd.read_csv(file))
    createPts(ptDict)

def write():
    outputfile = "Output/aki.csv"
    getNumAKI().to_csv(outputfile, index=False)

if __name__ == "__main__":
    files = glob.glob("Input/*.csv")
    print("The following files are proccessed: ", files)
    read(files)
    print("Processed", getNumPatients(), "patients, and", getNumCrs(), "laboratory values.")
    print("Writing results to disk")    
    write()
    input("Done, press enter to exit...")
