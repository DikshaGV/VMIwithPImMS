#TO GET CORRECTED DATA FOR PImMS VMI ANALSYIS FOR INDIVIDUAL PImMS RUN NUMBERS 
#Two functions are there 1) FEL is involved         2) FEL is not involved
#Input needed: DAQ run numbers (corresponding to the PImMS run numbers), HDF5 files path, `input.inp' file to add paths of the data required  
#Output: `alldata_<run number>rn.txt', where <run number> is the DAQ RUN NUMBER



#importing API
import runDataLoader as rdl
import beamtimeDataSettingsManager as bdsm
import fancyAnalysis as fa
#importing other useful packages
import numpy as np
import re
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy.interpolate import interp1d
from datetime import datetime 



def FEL_data(RunNumbers):
    for file in RunNumbers:
        dset = bdsm.SettingsForDataExtraction()
        dset.addMultipleChannelsFromConfigFile('input.inp')
        dataset = rdl.getRawCleanData(file,"/asap3/flash/gpfs/bl1/2018/data/11003927/raw/hdf/online-2",dset)
        dataset=fa.doTheJitterCorrection(dataset, PPDelay="PPdelay", BAMcorr="BAMcorr", PumpMask = None,ProbeMask =None, DumpInformation=True, ApplyJitterShift=False)

        PumpMask="FELShutter_markup"
        ProbeMask="ProbeEnergy_markup"
        TOFMS="TOF"
        MSType = "pump-probe"
        Mask = fa.getMaskAccordingToMode(dataset, TOFMS, PumpMask, ProbeMask, MSType)
        f = open("alldata_"+str(file)+"rn.txt", "w")
        data = dataset["ShotID"][Mask]
        print(np.shape(dataset["ShotID"][Mask]))
        for i in range(len(data)):
            f.write('%d\t%f\t%f\t%f\n'%(dataset["ShotID"][Mask][i],dataset["PPdelay"][Mask][i],dataset["PumpEnergy"][Mask][i],dataset["ProbeEnergy"][Mask][i])) #BunchID, Pump probe delay, PumpEnergy, ProbeEnergy





def NOFEL_data(RunNumbers):
    for file in RunNumbers:
        dset = bdsm.SettingsForDataExtraction()
        dset.addMultipleChannelsFromConfigFile('input_noFEL.inp')
        dataset = rdl.getRawCleanData(file,"/asap3/flash/gpfs/bl1/2018/data/11003927/raw/hdf/online-2",dset)

        TOFMS="TOF"
        f = open("alldata_"+str(file)+"rn.txt", "w")
        data = dataset["ShotID"]
        print(np.shape(dataset["ShotID"]))
        for i in range(len(data)):
            f.write('%d\t%f\t%f\t%f\n'%(dataset["ShotID"][i],dataset["PPdelay"][i],dataset["PumpEnergy"][i],dataset["ProbeEnergy"][i])) #BunchID, Pump probe delay, PumpEnergy, ProbeEnergy FOR OFFSHIFT 400-800 data fluorene


if __name__ == "__main__":
    
    RunNumbers = [21599,21600,21607,21608,21609,21610,21612,21613,21629] #Experimental DAQ run numbers
    FEL_data(RunNumbers)
    #NOFEL_data(RunNumbers)
    
    exit(0)
