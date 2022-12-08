#This script does many jobs :D
#1. It reads the.bin file which contains the data from PImMS camera for ion image and the TOF values for various frames
#2. It then clears the noisy pixels
#3. It can then select the particular ion fragment window
#4. Then it can bin the ion image according to the jitter corrected delay and signal corrected for FEL fluctuation 
#5. The binned ion image is then abel inverted and the radial distribution is stored
#6. The radial distribution is finally plotted against the delay time to get radial heatmaps

#READ THIS BEFORE RUNNING THIS SCRIPT
#You might want to change the run numbers you want to work for in the script main function.
#You might want to change the dimentions of the TOF array, for now value till 5000 are used because there is nothing in the spectrum beyond 2500
#Choose functions you want to work with, wisely!!
#To change the fragmnt ion,change the time window and the center of ion image

import numpy as np
import matplotlib.pyplot as plt
import abel
from abel import transform
import sys
import collections
from abel import linbasex
from mpl_toolkits import mplot3d
from scipy.interpolate import griddata
from datetime import datetime
from matplotlib import ticker
from matplotlib import cm
from PyInquirer import style_from_dict, Token, prompt, Separator
from types import SimpleNamespace


#---------------------------
#Load file NO FEL
#---------------------------
def read_bin_file_noFEL(file): #File: binary file, c,d: TOF gating, n: no. of photons (by default 1)          

    i = 0
    frames = 0
    ion = np.zeros(shape=(325, 325))
    tof = np.zeros(shape=(5000))

    path = "/asap3/flash/gpfs/bl1/2018/data/11003927/raw/pimms/"
    file1 = path + str(file)+".bin"  #Needed for the ion images
    file2 = path + str(file)+"-settings.txt" #Needed for ppdelay values. Required only for OFF-SHIFT data
    f = open(file1, "rb")
    ppdelay_file = np.loadtxt(file2)
    ppdelay_count=np.shape(ppdelay_file)

    ppdelay=[]
    ion_matrix = []
    tof_matrix = []

    while True and frames<50: 
    #while True:  
        f.seek(i, 0)
        params = np.fromfile(f, dtype=np.int32, count=2)
        if (np.shape(params) == (0,)):
            break
        data = np.fromfile(f, dtype=np.uint16, count=params[0] * params[1])
        i = f.tell()
        m = int(len(data) / 3)
        data_new = np.reshape(data, (m, 3))

        ion = np.zeros(shape=(325, 325))
        tof = np.zeros(shape=(5000))       
        for k in range(m):
            if data_new[k,2] > c and data_new[k,2] < d :          #Time window for fragment               
                ion[data_new[k, 0]][data_new[k, 1]] +=1 
                tof[data_new[k, 2]] += 1
        ion_matrix.append(ion)
        tof_matrix.append(tof)
        ppdelay.append(ppdelay_file[frames,3])

        print(frames)
        frames = frames + 1
    ion_matrix = np.array(ion_matrix)
    tof_matrix = np.array(tof_matrix)
    ppdelay = np.array(ppdelay)


    return(ion_matrix,tof_matrix,ppdelay)

#---------------------------
#Load file FEL
#---------------------------
def read_bin_file(file): #File: binary file, c,d: TOF gating, n: no. of photons (by default 1)          

    rem_frame=np.loadtxt("frames_to_remove_"+str(file)+".txt",dtype=int) #this file stores the name of the frame that is to be removed
    values=np.loadtxt("frame_BID_FELpower_corrdelay_"+str(file)+".txt") #This file store in each row "Frame number, BunchId, FELpower, delay value
    mean_FELpow = np.mean(values[:,2])
    #print(mean_FELpow)
    i = 0
    frame = 0
    ion = np.zeros(shape=(325, 325))
    tof = np.zeros(shape=(5000))
    path = "/asap3/flash/gpfs/bl1/2018/data/11003927/raw/pimms/"
    file1 = path + str(file)+".bin"
    f = open(file1, "rb")
    total_frames = 0

    FELpower=[]
    ion_int=[]
    FRAME=[]
    lst=[]
    delay=[]
    ion_matrix = []
    tof_matrix = []

    while True and total_frames<30: 
    #while True:  
        f.seek(i, 0)
        params = np.fromfile(f, dtype=np.int32, count=2)
        if (np.shape(params) == (0,)):
            break
        data = np.fromfile(f, dtype=np.uint16, count=params[0] * params[1])
        i = f.tell()
        if frame not in rem_frame:
            m = int(len(data) / 3)
            data_new = np.reshape(data, (m, 3))

            ion = np.zeros(shape=(325, 325))
            tof = np.zeros(shape=(5000))       
            for k in range(m): 
                if data_new[k,2] > c and data_new[k,2] < d :          #Time window for fragment
                    ion[data_new[k, 0]][data_new[k, 1]] +=1 
                    tof[data_new[k, 2]] += 1

            x = mean_FELpow/np.power(values[total_frames,2],n) #correction for FEL fluctuation
            ion_matrix.append(ion*x)
            tof_matrix.append(tof*x)

            delay.append(values[total_frames,3]) #delay already corrected for jitter
            FRAME.append(total_frames)
            total_frames=total_frames+1


        print(total_frames)
        frame = frame + 1
    ion_matrix = np.array(ion_matrix)
    tof_matrix = np.array(tof_matrix)
    ppdelay = np.array(delay)
    return(ion_matrix,tof_matrix,ppdelay)
    
#---------------------------
#get the tOF anf VMI of one ion or many
#---------------------------    
def getTOF_VMI():

    ion_matrix_final = []
    tof_final = []
    ppdelay_final=[]

    for i in R:
        if (FEL==True):
            ion_matrix,tof,ppdelay=read_bin_file(i)
        else:
            ion_matrix,tof,ppdelay=read_bin_file_noFEL(i)
        ion_matrix_final.append(ion_matrix)
        tof_final.append(tof)
        ppdelay_final.append(ppdelay)

    ion_matrix=np.vstack(ion_matrix_final)
    tof=np.vstack(tof_final)
    ppdelay=np.hstack(ppdelay_final)
    print("Shape of ion matrix "+str(np.shape(ion_matrix)))
    print("Shape of TOF array "+str(np.shape(tof))) 
    print("Shape of PPdelay array "+str(np.shape(ppdelay))) 
    print(ppdelay)

    #SAVING .dat files
    tof_sum=sum(tof) 
    np.savetxt("TOF_"+str(R)+".txt",tof_sum)

    #Plotting the preliminary VMI image
    ion_matrix_sum=sum(ion_matrix)
    np.savetxt("raw_Image_"+str(dict)+"_"+str(R)+"_.txt",ion_matrix_sum)
    fig, axs = plt.subplots()
    axs.set_xlabel('x-position')
    axs.set_ylabel('y-position')
    axs.matshow(ion_matrix_sum,origin='lower', extent=(0, 324, 0, 324))
    plt.show()
    return(ion_matrix,tof,ppdelay,ion_matrix_sum)

#---------------------------
#gfind the center
#--------------------------- 
def findcenter(ion_matrix_sum):
    questions = [
        {
            'type': 'input',
            'name': 'center_x',
            'message': 'Center of ion image (x value)',
            'default': '136.9'
        },
        {
            'type': 'input',
            'name': 'center_y',
            'message': 'Center of ion image (y value)',
            'default': '158.0'
        },
        {
            'type': 'input',
            'name': 'radius',
            'message': 'radius of the circle',
            'default': '65.0'
        }]
    answers = prompt(questions)
    m = SimpleNamespace(**answers)
    center=(float(m.center_x),float(m.center_y))
    radius = float(m.radius)
    circle = plt.Circle(center,radius,fill=False)
    fig, axs = plt.subplots()
    axs.matshow(ion_matrix_sum,origin='lower')
    axs.add_patch(circle)
    axs.set_xlabel('x-position')
    axs.set_ylabel('y-position')
    axs.set_title('Ion_Image')
    plt.show(block=False)
    print("Press any keyboard key to continue!")
    while True:
        if plt.waitforbuttonpress():
            break

    questions = [
        {
            'type': 'input',
            'name': 'choice',
            'message': '1. Satisfied? \n 2. Not satisfied\n 3. Exit and return th default values\nEnter your choice (1/2/3) : ',
            'default': '2.'
        }]
    answers = prompt(questions)
    m = SimpleNamespace(**answers)
    flag = m.choice
    if flag in {'1','1.'}:
        return(center)
    elif flag in {'3','3.'}:
        return((136.9,158.0))
    else:
        findcenter(ion_matrix_sum) 
    
    return(center)

#---------------------------
#save the centered image in a text file
#---------------------------     
def getcenteredVMI_static(ion_matrix_sum):
    center=addcenter(ion_matrix_sum_dummy) 
    center_new=tuple(reversed(center))
    cIM = abel.tools.center.center_image(ion_matrix_sum_dummy, method='com', odd_size=True, square=False, axes=(0, 1), order=3, verbose=False, center=center_new) #finding the centered image
    np.savetxt("centered_Image_"+str(dict)+"_"+str(R)+"_.txt",cIM)
    return(center)
    
#---------------------------
#shows the raw, abel inveted and radial intergration data, save the Abel inverted static image in a text file with its radial distribution. 
#--------------------------- 
def abel_invert(center,ion_matrix_sum):
    center_new=tuple(reversed(center))
    res = abel.transform.Transform(ion_matrix_sum, direction='inverse',center = center_new,recast_as_float64=True, symmetry_axis = (0, 1),method='onion_peeling',angular_integration=True)
    dist=res.angular_integration
    inverse_abel=res.transform

    fig, axs = plt.subplots(1, 3)

    pos=axs[0].imshow(ion_matrix_sum,origin='lower', extent=(0, 324, 0, 324))
    axs[0].set_xlabel('x-position')
    axs[0].set_ylabel('y-position')
    axs[0].set_title('Ion image')
    fig.colorbar(pos, ax=axs[0],shrink=0.5)
    pos2=axs[1].imshow(inverse_abel,origin='lower', cmap='viridis', extent=(0, 324, 0, 324))
    axs[1].set_xlabel('x-position')
    axs[1].set_ylabel('y-position')
    axs[1].set_title('Inverse Abel Transform\n')
    fig.colorbar(pos2, ax=axs[1],shrink=0.5)
    dist=np.array(dist)
    axs[2].plot(dist[0,:],dist[1,:])
    plt.show(block=False)
    print("Press any keyboard key to continue!")
    while True:
        if plt.waitforbuttonpress():
            break  
    #Saving the radial distribution in text file and Abel inverted ion image in .svg format
    dist2=np.transpose(dist) 
    np.savetxt("radial_dist_averageddelay_"+str(dict)+"_"+str(R)+".dat",dist2)
    np.savetxt("Abel_inverted_Image_sum"+str(dict)+"_"+str(R)+"_.txt",inverse_abel)
    return(dist2)

#---------------------------
#Plot the iamges.
#--------------------------- 

def plot_VMIprofile(dist_static,pix_to_mom):
    fig = plt.figure(figsize=(8,5))
    dist_static[:,0]=(dist_static[:,0]*pix_to_mom)
    plt.plot(dist_static[:,0],dist_static[:,1],linewidth=1.5,label = dict)
    plt.show()

#---------------------------
#Plot the iamges.
#---------------------------             
def  plot_images(file_name):
    
    IM = np.loadtxt(file_name)
    #For normalizing the image
    #p = np.max(IM)
    #print(p)
    #IM = IM/p
    
    fig, axs = plt.subplots()
    pos=plt.imshow(IM,origin='lower', cmap='gist_stern',extent=(0, 324, 0, 324))
    fig.colorbar(pos)
    plt.xticks(fontsize = '18')
    plt.yticks(fontsize = '18')
    plt.xlabel("x-position",fontsize='18')
    plt.ylabel("y-position",fontsize='18')
    #plt.title("Abel_inverted_image",fontsize='20')
    #plt.legend(loc='upper left',fontsize = '12')
    plt.savefig(str(dict)+".svg", dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=True, bbox_inches='tight', pad_inches=0.1,
            frameon=None, metadata=None)
    plt.show(block=False)
    print("Press any keyboard key to continue!")
    while True:
        if plt.waitforbuttonpress():
            break
    
#---------------------------
#Binning of the data based on ppdelay
#--------------------------- 
def bin_data(ion_matrix,ppdelay,bin_value):
    ini = np.min(ppdelay)
    fin = np.max(ppdelay)
    bin = bin_value

    ion_matrix_bin=[]
    ppdelay_bin=[]
    j = ini
    while j <= fin :
        ion_dummy = np.zeros(shape=(325, 325)) 
        points = 0
        for i in range(len(ppdelay)):
            if ppdelay[i]>=j and ppdelay[i]<j+bin:
                ion_dummy += ion_matrix[i]                
                #ppdelay_dummy += ppdelay[i]
                points +=1
        #ion_matrix_bin.append(ion_dummy)
        #ppdelay_bin.append(j+(bin/2.0))
        if points>0 :
            ion_matrix_bin.append(ion_dummy/points)
            ppdelay_bin.append(j+(bin/2.0))
        #    ion_matrix_bin.append(ion_dummy/points)
        #    ppdelay_bin.append(ppdelay_dummy/points)
        j=j+bin
   
    ion_matrix_bin=np.array(ion_matrix_bin)
    ppdelay_bin=np.array(ppdelay_bin)
    np.savetxt("ppdelay_bin_stack_"+str(dict)+"_"+str(R)+"_.txt",ppdelay_bin)
    return(ion_matrix_bin,ppdelay_bin)

#--------------------------------------------
#Abel inversion of binned matrix gives the radial distribution binned 
#----------------------------------------------------------------
def abel_inversion_bin_data(ion_matrix_bin,center):
    center_new=tuple(reversed(center))

    inverse_abel = []
    dist = []
    for r in range(len(ion_matrix_bin)):
        #print(":")
        print(r)
        res = abel.transform.Transform(ion_matrix_bin[r], direction='inverse',center = center_new, symmetry_axis = (0, 1),recast_as_float64=True, method='onion_peeling',angular_integration=True)
        #res = abel.transform.Transform(ion_matrix_bin[r], direction='inverse',center = center_new,recast_as_float64=True, method='onion_peeling',angular_integration=True)
        inverse_abel.append(res.transform)
        dist.append(res.angular_integration)
    
    inverse_abel=np.array(inverse_abel)
    dist=np.array(dist)
    dist_bin_stack=np.vstack(dist)
    np.savetxt("dist_bin_stack_"+str(dict)+"_"+str(R)+"_.txt",dist_bin_stack)
    
    #To check the data after binning and abel inverted here  
    fig, axs = plt.subplots()
    ion_matrix_sum=sum(inverse_abel)
    axs.imshow(ion_matrix_sum,origin='lower', extent=(0, 324, 0, 324))
    axs.set_xlabel('x-position')
    axs.set_ylabel('y-position')
    axs.set_title('Sum of Abel invertd images after binning')
    plt.show(block=False)
    print("Press any keyboard key to continue!")
    while True:
        if plt.waitforbuttonpress():
            break
    return(dist)

#--------------------------------------------
#Arranging the data dist_bin and ppdelay_bin to create heatmap in origin
#--------------------------------------------
def arrange(ppdelay_bin,dist_bin,pix_to_mom): 
    ppdelay_len=len(ppdelay_bin)
    dist_len=len(dist_bin[0,0,:])
    dist_time=[]
    for r in range(len(ppdelay_bin)):
        a=[]
        for i in range(dist_len):
            a.append(ppdelay_bin[r])
        a=np.transpose(np.array(a))
        dist_time.append(np.transpose((np.vstack((dist_bin[r],a)))))    
    dist_time=np.array(dist_time)
    
    #dist_time=np.transpose(dist_time)
    print(np.shape(dist_time))

    dist_time_final=dist_time[0]   
    for r in range(ppdelay_len-1):
        dist_time_final=np.vstack((dist_time_final,dist_time[r+1]))
    dist_time_final[:, [0,1]] = dist_time_final[:, [1,0]] 
    dist_time_final[:, [0,2]] = dist_time_final[:, [2,0]]
    print("finally the data is ready for ORIGIN to plot given by file \"time_pixel_intensity_....dat\"")
    print(np.shape(dist_time_final))
    np.savetxt("delaytime_pixel_intensity_"+str(dict)+"_"+str(R)+"_.dat",dist_time_final)
    
    #Momentum calibration
    f = open("delaytime_momentum_intensity_"+str(dict)+"_"+str(bin_value)+"_"+str(R)+"_.dat", "w")
    for i in range(len(dist_time_final)):
        f.write('%r\t%r\t%r\n'%(dist_time_final[i,0],(dist_time_final[i,1]*pix_to_mom)/137.0,dist_time_final[i,2])) #Delay momentum Intensity
    f.close()

    return(dist_time_final)
#-----------------------------------------
def plot_ionyieldvsppdelay(pix_int,ppdelay_bin,ion,px0,px1,px2,px3):  
    dist=pix_int
    ini,fin,a,b=px0,px3,px1,px2
    # For integration of  {x = dist[i,0,:] (pixels)} and { y = dist[i,1,:] (intensity) } where i is the ppdelay ==> len(dist) = ppdelay length
    # The number of points to integrate is length of x or y
    
    total_points = len(dist[0,0,:])
    integral=[]
    integral_low,integral_med,integral_high=[],[],[]
    for j in range(len(dist)):
        sum,sum_low,low_size,sum_med,med_size,sum_high,high_size=0.0,0.0,0.0,0.0,0.0,0.0,0.0
        for i in range(total_points):
            sum+=dist[j,1,i]
            if dist[j,0,i] < a and dist[j,0,i]>ini:
                sum_low+=dist[j,1,i]             #To find the sum in the bins
                low_size+=1                      #no. of points those are summed
            elif dist[j,0,i] <b and dist[j,0,i]>a:
                sum_med+=dist[j,1,i]
                med_size+=1
            elif dist[j,0,i]>b and dist[j,0,i]<fin:
                sum_high+=dist[j,1,i]
                high_size+=1
        all_points = low_size+med_size+high_size
        sum_low=sum_low/all_points               # finding the average by dividing the sum of points to the no. of the points
        sum_med=sum_med/all_points
        sum_high=sum_high/all_points
        sum = sum/all_points
       
        integral.append(sum)
        integral_low.append(sum_low)
        integral_med.append(sum_med)
        integral_high.append(sum_high)
    integral=np.array(integral)
    integral_low=np.array(integral_low)
    integral_med=np.array(integral_med)
    integral_high=np.array(integral_high)
    
    f =open(str(ion)+"_"+str(R)+"_ppiy_new.dat","w")
    for i in range(len(ppdelay_bin)):
        f.write('%f\t%f\t%f\t%f\t%f\n'%(ppdelay_bin[i],integral_low[i],integral_med[i],integral_high[i],integral[i]))
    f.close()

    #plotting the curves
    fig,axs = plt.subplots(1,3)
    axs[0].plot(ppdelay_bin,integral_low)
    axs[0].set_title(ion+"(1,0)")
    axs[1].plot(ppdelay_bin,integral_med)
    axs[1].set_title(ion+"(1,1)")
    axs[2].plot(ppdelay_bin,integral_high)
    axs[2].set_title(ion+"(1,2)")
    fig,ax = plt.subplots()
    ax.plot(ppdelay_bin,integral)
    plt.show()

 
def IonyieldvsPumpprobe(ions,minimum,maximum,band_min,band_max):
    for i in range (len(ions)):
        ion,px0,px1,px2,px3 = ions[i],minimum[i],band_min[i],band_max[i],maximum[i]
        pix_int = np.loadtxt("dist_bin_stack_"+str(ion)+"_"+str(R)+"_.txt")    
        ppdelay_bin = np.loadtxt("ppdelay_bin_stack_"+str(dict)+"_"+str(R)+"_.txt")
        l = len(pix_int)
        pix_int = np.vsplit(pix_int,l/2)
        pix_int = np.array(pix_int) 
        print(np.shape(pix_int))
        #print(pix_int[0])
        plot_ionyieldvsppdelay(pix_int,ppdelay_bin,ion,px0,px1,px2,px3)
        #print (values[0])
        
        
        
if __name__=="__main__":


    #careful with saving the images, takes alot of memory :-)

    global R,dict,c,d,n,FEL,bin_value,pix_to_mom #R:PImMS run number
                            #dict: name of the fragment ion
 		                    #[c,d]: gating for the ion based on TOF
                            #n: no. of photons involve
                            #FEL: flag  (True: FEL used, False: FEL not used)
                            #bin_value: binning on ppdely time
                            #pix_to_mom: pixel to momentum calibration factor for the ions                             

    FEL=True 
    R = [275,276]
    bin_value=0.05
    c,d,dict,n,center,pix_to_mom=2116,2125,"something",1,(162,162),135/137


    #######################
    #   LOADING DATA      #
    ################################################################################################################
    # 1. getTOF_VMI() loads and returns the raw ion_matrix, tof, ppdelay arrays                                   
    #    Input: uses only the global paramters (Run number, fragment ion: name and tof window, number of photons) 
    #    outputs: TOF and VMI averaged over all time delays (also, can be called a part of static analysis)       
    ################################################################################################################

    ion_matrix,tof,ppdelay,ion_matrix_sum=getTOF_VMI()
    
    ##################################################################################################
    #   STATIC analysis (ion_matrix_sum can be read from a file to save time from loading the data)  #
    ##################################################################################################
    # 1. Find the center of the image "THIS FUNCTION IS REDUNDANT IF YOU KNOW THE CENTER"                
    #    Input: raw VMI file and global paramters 
    #    Outputs: center values in the (x,y) tupple fromat, not in (row,columns) tuple format so need to be reversed  
    #             center VMI image data file  
    #
    # 2. Abel_invert the delay averaged
    #    Input: raw VMI file,center, and global paramters
    #    Outputs: Abel inverted image and its radial distribution in a text file  
    #    Plots: Raw image, Abel inverted image, radial distribution. 
    # 3. Plotting momentum profile of the fragment ions
    #    Plots and save the any Image in matrix form, can tune the color sclae etc, plotting format etc.      
    # 4. Plotting function  OPTIONAL TO USE
    #    Plots any Image in matrix form, can tune the color sclae etc, plotting format etc. 
    ###########################################################################################################
     
    #1.
    #ion_matrix_sum=np.loadtxt("raw_"+str(dict)+"_"+str(R)+"_.txt") #or you can pass the ion_matrix_sum_directly that you get from the getTOF_VMI()
    #center=getcenteredVMI_static(ion_matrix_sum) #Find the center of the image "THIS FUNCTION IS REDUNDANT IF YOU KNOW THE CENTER"
    #2.
    #dist_static=abel_invert(center,ion_matrix_sum) 
    #3.
    #dist_static=np.loadtxt("radial_dist_averageddelay_"+str(dict)+"_"+str(R)+".dat")
    #plot_VMIprofile(dist_static,pix_to_mom)
    #4.
    #file_name="raw_"+str(dict)+"_"+str(R)+"_.txt"
    #plot_images(file_name)
    
    ##################################################################################################
    #   DYNAMIC analysis                                                                             #
    ##################################################################################################
    # 1. Bin the data                
    #    Input: raw ion_matrix list per BID, ppdely values per BID, bin_value needed, and global paramters 
    #    Outputs: binned ion_matrix, and ppdelay_bin values   
    #
    # 2. Abel_invert the the binned ion images and calulate the radial distribution.
    #    Input: just the binned ion matrix, center, and global paramters
    #    Outputs: radial distribution binned (also saved as tet file) 
    #    Plots: [CHECK]sum of abel inverted binned images to check whether the binning ruined the data... It should look exactly the same as the static Abel inverted image. !!!!!    
    # 3. Arranging the file such that 2 separate data ppdelay_bin and dist_bin (pixel vs intensity in dimension (2,230)) per delayed dimension for example: (binneddelay_array_size,2,230) 
    #    Input: binned radial distribution and the binned ppdelay.
    #    Output: data in the 3 column as ppdelay,pixel,intensity and also as ppdelay,momentum,intensity
    ###########################################################################################################
    
    ion_matrix_bin,ppdelay_bin=bin_data(ion_matrix,ppdelay,bin_value)
    dist_bin=abel_inversion_bin_data(ion_matrix_bin,center)
    arrange(ppdelay_bin,dist_bin,pix_to_mom)
    
    ##################################################################################################
    #   POST-PROCESSING functions                                                                             #
    ##################################################################################################
    # 1. Ionyieldvspumpprobe delay                
    #    Input: distribtuion_binned, pixel values of the fragments for different momentum regions.
    #    Outputs: ion yield v/s pump probe delay for low momentum, medium momentum, high momentum, all the momentum regions.   
    ###########################################################################################################
    #ions = ["C5hx++", "C7hx++", "C11hx++"]
    #minimum,maximum = [3.5,3.5,3.5],[104,94,55] #new values
    #band_min = [42.3,38,18.5]
    #band_max = [100,64.2,37.5]
    ions = [dict]
    minimum,maximum = [3.5],[104] #new values
    band_min = [42.3]
    band_max = [100]
    IonyieldvsPumpprobe(ions,minimum,maximum,band_min,band_max)
  