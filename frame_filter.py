#This file does the syncing between the PImMS data and the HDF data suing common Bunch Id as tag.
#Input files: "alldata_<rn>rn.txt" and "-settings-corrected-v12.txt"
#Output files: "frame_BID_FELpower_corrdelay_<rn>.txt" and "frames_to_remove_<rn>.txt"

import numpy as np


#each frame has one BunchId

mcp = [21599,21600,21607,21608,21609,21610,21612,21613,21629]
pimms = [275,276,283,284,285,286,288,289,304]
for h in range(len(pimms)):
    file,daqfile = pimms[h],mcp[h]
    path2 = "/asap3/flash/gpfs/bl1/2018/data/11003927/processed/pimms/"
    file2 = path2 + str(file)+"-settings-corrected-v12.txt"

#INPUT FILES    
    pimms_data=np.loadtxt(file2)
    daq_data=np.loadtxt("alldata_"+str(daqfile)+"rn.txt")
#OUTPUT FILES
    f1=open("frame_BID_FELpower_corrdelay_"+str(file)+".txt","w")
    f2=open("frames_to_remove_"+str(file)+".txt","w")

#print(pimms_data[1,10])
#print(daq_data[:,0])

    for i in range(len(pimms_data)):
    
        for j in range(len(daq_data)):
        
            if (pimms_data[i,10] == daq_data[j,0]) and pimms_data[i,11] != 0.0:
                print(i,file)
                f1.write('%d\t%d\t%f\t%f\n'%(i,pimms_data[i,10],pimms_data[i,11],daq_data[j,1]))
                k=-1
                break
            else:
                k=i
        if (k!=-1):
            f2.write('%i\n'%k)