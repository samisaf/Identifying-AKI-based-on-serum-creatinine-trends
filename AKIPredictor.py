# -*- coding: utf-8 -*-
"""@author: Sami Safadi"""

from patient import Patient
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob

demDict = pd.DataFrame()
ptDict = {} # Dictionary containing all lab values
Patients = {} # Dictonary, keys are MRN, values are Patient objects
        
def getNumPatients():
    return len(Patients)
    
def getNumCrs():
    return sum([Patients[p].crs.size for p in Patients])
     
def createPtDict(df: pd.DataFrame):
    global ptDict 
    uniqueMRN = np.unique(df.ix[:, 0])
    for i in uniqueMRN: 
        if not(i in ptDict.keys()): ptDict[i] = dict()
    for i in range(len(df)):
        mrn = df.ix[i, 0]
        creatinine = df.ix[i, 1]
        date = df.ix[i, 2]
        if not(np.isnan(mrn) or np.isnan(creatinine)): 
            ptDict[int(mrn)][str(date)] = float(creatinine)
    return ptDict

def createPts(patients: dict):
    global Patients
    for key in patients:
        crs = patients[key]
        mrn = key
        age, gender, race = 0, 0, 0
        if any(demDict.MRN == 20):
            index = demDict.index[demDict.MRN == key][0]
            age = demDict.Age[index]
            gender = demDict.Gender[index]
            race = demDict.Race[index]
        if len(crs) > 0: Patients[mrn] = Patient(mrn, crs, age, gender, race)

def getPlots(keys: list, savetofile = False):
    global Patients
    for key in keys:
        p = Patients[key]
        plt.figure()
        p.plot()
        if savetofile: 
            outputfile = "Output/Graphs/" + str(key) + ".png"
            plt.savefig(outputfile)
        else: plt.show()

def getAKIDates(keys: list, savetofile = False):
    global Patients
    for key in keys:
        akiDates = Patients[key].aki
        if savetofile:
            outputfile = "Output/Dates/" + str(key) + ".csv"
            if akiDates.size > 0: akiDates.to_csv(outputfile)
        else:
            print(akiDates)

def getNumAKI(keys: list, savetofile = False):
    global Patients
    mrn, aki, baseCr, egfr = [], [], [], []
    for key in keys:
        mrn.append(Patients[key].mrn)
        aki.append(Patients[key].aki.size)
        baseCr.append(Patients[key].baseCr)
        egfr.append(Patients[key].egfr)
    di = {'MRN': mrn, 'baseCr': baseCr, 'numAKI': aki, 'eGFR' : egfr}
    df = pd.DataFrame(di)
    df['anyAKI'] = df.numAKI > 0
    outputfile = "Output/aki.csv"
    if savetofile: df.to_csv(outputfile, index = False)
    else: print(df)

def readDemographics(file: str):
    global demDict
    demDict = pd.read_csv(file)
        
def readLabs(files:[str]):
    for file in files: 
        createPtDict(pd.read_csv(file))
    createPts(ptDict)
    
def write():
    global Patients
    getNumAKI(Patients.keys(), savetofile = True)
    getAKIDates(Patients.keys(), savetofile = True)
    getPlots(Patients.keys(), savetofile = True)
    
if __name__ == "__main__":
    demographicsFile = "Input/Demographics.csv"
    readDemographics(demographicsFile)
    
    labFiles = glob.glob("Input/Labs*.csv")
    print("The following files are proccessed: ", labFiles)
    readLabs(labFiles)
    
    message = "Processed {} patients, and {} laboratory values"\
        .format(getNumPatients(), getNumCrs(),)
    print(message)
    print("Writing results to disk")    
    write()
    input("Done, press enter to exit...")
