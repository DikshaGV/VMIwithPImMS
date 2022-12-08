#! /usr/bin/env python

# importing libraries for FLASH analysis
import runDataLoader as rdl
import beamtimeDataSettingsManager as bdsm
import fancyAnalysis as fa

# importing libraries for plotting
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d

#importing libraries for building command line interfaces
import argparse
from PyInquirer import style_from_dict, Token, prompt, Separator
from examples import custom_style_2
from types import SimpleNamespace 

matplotlib.rc('font', size='20')

   
#************************************************************************************************************************
#Step 3
def MenuAnalysis():
    
    #************************************************************
    #predefined functions need to use the MenuAnalysis() function
    def twoDMS():
        ppd,ms = fa.get2dMassSpectrum(res, PumpMask="FELShutter_markup", ProbeMask="ProbeEnergy_markup", TOFMS="TOF", PPDelay="PPdelay", NumOfBins = 20, OutputForGnuplot="2d_tofms_power.tmp",
                             PumpEnergy="PumpEnergy",  PumpNumOfPhotons = 1, 
                             ProbeEnergy=None, ProbeNumOfPhotons = 1, RemoveBackground=True)


        ma, mdms = fa.convertMSFromTOFDomainToMByQDomain(ms,TOFCal,dx=0.05, ApplyJacobian=True, MaxMass=None, ExtrapolationType='cubic')
        
        outf=open('2d_mq.tmp', 'w')
        for i,t in enumerate(ppd):
            for j,m in enumerate(ma): 
                outf.write(" %15.10f %15.10f  %15.10f\n" % (t,m,mdms[i][j]))
            outf.write("\n")
        outf.close()

        #2D mass spectrum plot
        fig, ax = plt.subplots()
        ax = plt.axes(projection ='3d')
        x,y = np.meshgrid(ma,ppd)
        print(x.shape)
        print(y.shape)
        surf = ax.plot_surface(x,y,mdms, alpha=1,rstride=1,cstride=1,cmap=cm.coolwarm)
        ax.set(xlabel='\n\nM/Q, arb.units', ylabel='\n\npp delay', zlabel = '\n\nIntensity', title='2 D TOF Mass Spectrum for '+str(MStype)+' data\n\n')
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show(block=False)
        print("Press any keyboard key to continue!")
        while True:
            if plt.waitforbuttonpress():
                break

    def ppdelaycurves(fragment_window):
        data=fa.doTheJitterCorrection(res, PPDelay="PPdelay", BAMcorr="BAMcorr", PumpMask = None,ProbeMask =None, DumpInformation=True, ApplyJitterShift=False)

        fa.getIonYieldVsPumpProbeDelay(Data=data,
                                      Calibration=TOFCal, IonDefAndStore=fragment_window,
                                      PumpMask="FELShutter_markup", ProbeMask="ProbeEnergy_markup",
                                      PumpEnergy=None,  PumpNumOfPhotons = 1, 
                                      ProbeEnergy=None, ProbeNumOfPhotons = 1,
                                      TOFMS="TOF", PPDelay="PPdelay",
                                      NumOfBins=50, NormByJacobian=True, MinPtsNumber=100)
       
        fragment_window.saveIonYieldWindowToFile(fragment_window._name+"_"+ str(RunNumber) + "_window.dat")
        fragment_window.saveIonPPYieldToFile(fragment_window._name+"_"+ str(RunNumber) + "_ppiy.dat")
        ppyield = np.loadtxt(fragment_window._name +"_"+ str(RunNumber) +"_ppiy.dat",dtype=float)
        window = np.loadtxt(fragment_window._name +"_"+ str(RunNumber)+ "_window.dat",dtype=float)

        #plotting ppdelay curves
        fig, ax = plt.subplots(1,2)
        ax[0].errorbar(ppyield[:,0], ppyield[:,1], yerr=ppyield[:,3], xerr=ppyield[:,2])
        ax[0].set(xlabel='Pump-probe delay', ylabel='Ion Yield', title="pump-probe delay v/s Ion Yield curve ("+ fragment_window._name + ")\n")
        ax[0].grid()
        ax[1].plot(window[:,0],window[:,1])
        ax[1].set(xlabel='M/Q, arb.units', ylabel='Intensity', title="fragment window ("+ fragment_window._name + ")\n")
        plt.show(block=False)
        print("Press any keyboard key to continue!")
        while True:
            if plt.waitforbuttonpress():
                break


    def TofTofCovariance():
        res1 = fa.calcTOFvsTOFCovMtrx(Data=res, PumpMask="FELShutter_markup", ProbeMask="ProbeEnergy_markup",
                  TOFMS="TOF", PPDelay="PPdelay", PumpLaserIntensity="PumpEnergy", ProbeLaserIntensity=None,
                  PCovScale = 0.0,
                  t0 = None, tWidth = 0.0, BeforeT0 = True,
                  startTOF = 0, endTOF = None, startTOF2 = 0, endTOF2 = None,MSType=MStype)

        np.savetxt("pure_cov_full.dat", res1["covmtrx"])
        np.savetxt("corr_full.dat", res1["pumpcorr"])
  

        pure_cov=np.loadtxt("pure_cov_full.dat",dtype=float)
        pump_corr=np.loadtxt("corr_full.dat",dtype=float)
        return pure_cov,pump_corr


    def Plot_Cov(pure_cov,pump_corr,s=0.0,layout='CovCorr_1',cmap_style='viridis',c_min=None,c_max=None):
        corr_cov=pure_cov-(s*pump_corr) #corr_cov is the covariance matrix corrected for laser intensity fluctuation (stored in pump_corr)
        np.savetxt("corr_cov.dat", corr_cov)
        pcov_max,pcov_min = corr_cov.max(),corr_cov.min()

        if layout == 'CovCorr_1':
         
            #covariance plot layout type 1
            fig, axs = plt.subplots(1,2)          
            plot1=axs[0].imshow(corr_cov,origin='lower',cmap=cmap_style,vmin=c_min,vmax=c_max)
            plot2=axs[1].imshow(pump_corr,origin='lower',cmap=cmap_style,vmin=c_min,vmax=c_max)
            axs[0].set_title('Corrected TOF-TOF Covariance for s = '+str(s))
            axs[1].set_title('Correction matrix for laser intensity fluctuation')
            fig.colorbar(plot1, ax=axs[0],shrink=0.5)
            fig.colorbar(plot2, ax=axs[1],shrink=0.5)
            plt.show(block=False)
            print("Press any keyboard key to continue!")
            while True:
                if plt.waitforbuttonpress():
                    break


        elif layout == 'CovCorr_2':
            
            #covariance plot layout type 2
            fig, axs = plt.subplots()          
            plot1=axs.imshow(corr_cov,origin='lower',cmap=cmap_style,vmin=c_min,vmax=c_max)
            axs.set_title('Corrected TOF-TOF Covariance for s = '+str(s))
            fig, ax = plt.subplots() 
            plot2=ax.imshow(pump_corr,origin='lower',cmap=cmap_style,vmin=c_min,vmax=c_max)
            ax.set_title('Correction matrix for laser intensity fluctuation')
            fig.colorbar(plot1,ax=axs)
            fig.colorbar(plot2,ax=ax)
            plt.show(block=False)
            print("Press any keyboard key to continue!")
            while True:
                if plt.waitforbuttonpress():
                    break
        elif layout == 'CovOnly':
            fig, axs = plt.subplots()          
            plot=axs.imshow(corr_cov,origin='lower',cmap=cmap_style,vmin=c_min,vmax=c_max)
            axs.set_title('Corrected TOF-TOF Covariance for s = '+str(s))
            fig.colorbar(plot,ax=axs,shrink=1.0)
            plt.show(block=False)
            print("Press any keyboard key to continue!")
            while True:
                if plt.waitforbuttonpress():
                    break
        return pcov_max,pcov_min

    def customize(pure_cov,pump_corr,pcov_max,pcov_min):        
        flag='0'
        while flag not in {'N','n'}:       
            print("CUSTOMIZE YOUR TOF-TOF COVARIANCE RESULTS")
            questions = [

                {
                    'type': 'input',
                    'name': 'scaling_factor',
                    'message': '\nEnter the scaling value, s  (Note: s=0.0 for pure covariance and s!=0.0 for partial covariance) : ',
                    'default': '0.0'
                },
                {
                    'type': 'list',
                    'message': 'Select the layout of the figure',
                    'name': 'layout',
                    'choices': ['1. Two plots in one graphical window (Left panel : Covariance matrix; Right panel : Correction matrix',
                               '2. Two plots in two different graphical window : (First plot : Covariance matrix; Second plot : Correction matrix',
                               '3. Plot only covariance matrix'],
                },
                {
                    'type': 'input',
                    'name': 'cmap_style',
                    'message': '\nEnter the color scaling for the covariance plots \n'
                               '-Perceptually Uniform Sequential : [viridis, plasma, inferno, magma, cividis]\n'
                               '-Sequential : '
                               '[Greys, Purples, Blues, Greens, Oranges, Reds,'
                               'YlOrBr, YlOrRd, OrRd, PuRd, RdPu, BuPu,'
                               'GnBu, PuBu, YlGnBu, PuBuGn, BuGn, YlGn]\n'
                               'Note: For more options go to the link https://matplotlib.org/3.1.3/tutorials/colors/colormaps.html\n'
                               'Enter your choice : ',
                    'default': 'viridis'
                },
                {
                    'type': 'input',
                    'name': 'v_min',
                    'message': '\nEnter the value to represent minimum intensity color '
                               '(note: The minimum intensity value of covariance matrix is '+str(pcov_min)+') ): ',
                    'default': '0.0'
                },
                {
                    'type': 'input',
                    'name': 'v_max',
                    'message': '\nEnter the value to represent maximum intensity color '
                               '(note: The maximum intensity value of covariance matrix is '+str(pcov_max)+') ): ',
                    'default': '100.0'
                },]
             
            answers=prompt(questions,style=custom_style_2)
            s=float(answers['scaling_factor'])
            if answers['layout']=='1. Two plots in one graphical window (Left panel : Covariance matrix; Right panel : Correction matrix':
                layout = 'CovCorr_1' # 1 stands for plots in 1 figure 
            elif answers['layout']=='2. Two plots in two different graphical window : (First plot : Covariance matrix; Second plot : Correction matrix':
                layout = 'CovCorr_2' # 2 stands for plots in 2 figures
            elif answers['layout']=='3. Plot only covariance matrix':
                layout = 'CovOnly'
            cmap_style=answers['cmap_style']
            c_min=float(answers['v_min'])
            c_max=float(answers['v_max'])
            Plot_Cov(pure_cov,pump_corr,s,layout,cmap_style,c_min,c_max)
            
            print("Do you want to try different colors?(y/n)")
            flag=input()
            if flag in {'n','N'}:
                break

    
    flag = '0'
    while flag!='4.' or flag!='4':
        questions = [
            {
                'type': 'input',
                'name': 'save',
                'message': '\n1. 2 D Mass Spectrum \n2. TOF TOF Covariance\n3. Pump-probe delay curves\n4. Exit\nEnter your choice:  ',
                'default':''
            },]
        answers = prompt(questions, style=custom_style_2)
        n = SimpleNamespace(**answers)
        flag = n.save
        if flag in {'4.','4','exit'}:
            break

        if flag in {'1.','1','2 D Mass spectrum'}:
            twoDMS()

        elif flag in {'2.','2','TOF TOF covariance'}:
            pure_cov,pump_corr=TofTofCovariance()
            pcov_max,pcov_min=Plot_Cov(pure_cov,pump_corr,s=0.0,layout='CovCorr_1',cmap_style='viridis',c_min=None,c_max=None) 
            flag1='y'
            while flag1 in {'y','Y'}:
                questions = [

                    {
                        'type': 'input',
                        'name': 'choice',
                        'message': 'Do you want to customize your graphical TOF-TOF results?(y/n)\nEnter your choice : ',
                        'default': 'y'
                    }]
                answers = prompt(questions, style=custom_style_2) 
                n = SimpleNamespace(**answers)
                flag1=n.choice
                if flag1 in {'n','N'}: 
                    break  
                elif flag1 in {'Y','y'}:
                    customize(pure_cov,pump_corr,pcov_max,pcov_min)
                elif flag1 not in {'n','y','N','Y'}:
                    print("You entered wrong keyword. PLEASE ENTER (y/n)")
                
         
            

        elif flag in {'3.','3','Pump-probe delay curves'}:

            questions = [

                {
                    'type': 'input',
                    'name': 'fragment_name',
                    'message': '\nSpecify the ion name for pp delay curve in format \n',
                    'default': 'C2Hx'
                }, 
                {
                    'type': 'input',
                    'name': 'fragment_min',
                    'message': '\nSpecify the minimum mass in order to specify the window \n',
                    'default': '25.5'
                },
                {
                    'type': 'input',
                    'name': 'fragment_max',
                    'message': '\nSpecify the maximum mass in order to specify the window \n',
                    'default': '26.5'
                },]
            answers = prompt(questions, style=custom_style_2)
            n = SimpleNamespace(**answers)
            fragment_window=fa.IonicPumpProbeYield(n.fragment_name,float(n.fragment_min),float(n.fragment_max))
            ppdelaycurves(fragment_window)      


#*******************************************************************************************************************
#STEP 2   
#function to calibrate the time of flight to mass by charge ratio
def calib(ms):
    tof_spec=ms
    global TOFCal
    TOFCal = fa.TOFMSCalibration()


    #*************************************
    #pre defined functions need to be used by calib() 
    def plotting(maxis,ms):
        fig, ax = plt.subplots()
        ax.plot(maxis, ms)
        ax.set(xlabel='M/Q, arb.units', ylabel='Intensity',title='1 D TOF Mass Spectrum for '+str(MStype))
        ax.grid()
        plt.show(block=False)
        print("Press any keyboard key to continue!")
        while True:
            if plt.waitforbuttonpress():
                break
        
    
    def parceLine(line):
        words = line.split()
        if len(words)==3:
            return {"TOF": float(words[0]), "MbyQ": float(words[1]), "IonName": str(words[2])}
        elif len(words)==2:
            return {"TOF": float(words[0]), "MbyQ": float(words[1]), "IonName": None}
        else:
            print("Wrong format of the input")
            return False

    def save_calib(maxis,ms):
        questions = [
           {
                'type':'input',
                'name':'cal_save',
                'message':'save the current mass calibrated spectrum into text file? (y/n)',
                'default':'y'
            },]
        answers = prompt(questions, style=custom_style_2)
        n = SimpleNamespace(**answers)
        if n.cal_save in {'y','Y'}:
            np.savetxt(MolName+"_1D_tof_mass_calib_spect_"+str(MStype)+"_"+str(RunNumber)+".dat",np.stack((maxis,ms),axis=-1))
            print("Calibration saved into "+MolName+"_1D_mass_spectrum_"+str(MStype)+".dat successfully")
        elif n.cal_save in {'n','N'}:
            print("calibration not saved in the data file")
        else:
            print("wrong choice entered, calibration could not be saved into the text file!")


    def AddPoints():
        try:
            questions = [
                {
                    'type': 'input',
                    'name': 'line',
                    'message': 'Please give information of ion in the form [TOF] [M/Q] ([Ion Name]): ',
                    'default': ''
                },]
            answers = prompt(questions, style=custom_style_2)
            n = SimpleNamespace(**answers)
            ionInfo = parceLine(n.line)
            TOFCal.addIon(**ionInfo)
        except: 
            questions = [
                {
                    'type': 'input',
                    'name': 'line',
                    'message': 'Please give information of ion in the form [TOF] [M/Q] ([Ion Name]): ',
                    'default': ''
                },]
            answers = prompt(questions, style=custom_style_2)
            n = SimpleNamespace(**answers)
            ionInfo = parceLine(n.line)
            TOFCal.addIon(**ionInfo)

    #**********************************   
    #function calib starts here! 
    questions = [
        {
            'type': 'input',
            'name': 'choice',
            'message': 'Calibration (Please enter 1 or 2):\n1. Add two calibration points\n2. Use Calibration File for calibration\nEnter your choice: ',
            'default': '2'
        },]
    answers = prompt(questions, style=custom_style_2)
    m = SimpleNamespace(**answers)
    
    if m.choice in {'1.','1','Add calibration points'}:
        print("\n\nFor the First Ion")
        AddPoints()
        print("\n\nFor the Second Ion")
        AddPoints()

        TOFCal.display_ions()
 
        maxis, ms = fa.convertMSFromTOFDomainToMByQDomain(tof_spec, TOFCal, dx=0.05, ApplyJacobian=True, MaxMass=500., ExtrapolationType='cubic')

        #Mass spectrum plot
        plotting(maxis,ms)

        TOFCal.display_ions()
        save_calib(maxis,ms)

        flag = '0'
        while flag!='4.' or flag!='4':
            questions = [
                {
                    'type': 'input',
                    'name': 'save',
                    'message': '\n1. Remove ions from the library\n2. Add more ions?\n3. Save the existing ions in a txt file CalFile.inp!\n4. Save calibrated mass spectrum with current calibration points into text file \n5. Exit\nEnter your choice:  ',
                    'default': ''
                },]
            answers = prompt(questions, style=custom_style_2)
            n = SimpleNamespace(**answers)
            flag = n.save
            if flag in {'5.','5','exit'}:
                break
            if flag in {'1.','1','remove'}:
                questions = [
                    {
                        'type': 'input',
                        'name': 'IonName',
                        'message': 'enter the name of the Ion you want to remove ',
                        'default': ''
                    },]
                answers = prompt(questions, style=custom_style_2)
                n = SimpleNamespace(**answers)
                TOFCal.rmIon(n.IonName)
                TOFCal.display_ions()

            elif flag in {'2.','2','Add'}:
                AddPoints()
                TOFCal.display_ions()
                maxis, ms = fa.convertMSFromTOFDomainToMByQDomain(tof_spec, TOFCal, dx=0.05, ApplyJacobian=True, MaxMass=1000., ExtrapolationType='cubic')
                
                #Mass spectrum plot               
                plotting(maxis,ms)

            elif flag in {'3.','3','save'}:
                TOFCal.display_ions()
                TOFCal.saveCalibrationFile('CalFile.inp')

            elif flag in {'4.','4','save_calib'}:
                
                maxis, ms = fa.convertMSFromTOFDomainToMByQDomain(tof_spec, TOFCal, dx=0.05, ApplyJacobian=True, MaxMass=1000., ExtrapolationType='cubic')
                TOFCal.display_ions()
                save_calib(maxis,ms)
                                          

    elif m.choice in {'2.','2',' Use Calibration File for calibration'}:
        TOFCal.readCalibrationFile("CalFile.inp")
        maxis, ms = fa.convertMSFromTOFDomainToMByQDomain(tof_spec, TOFCal, dx=0.05, ApplyJacobian=True, MaxMass=1000., ExtrapolationType='cubic')
        
        #Mass spectrum plot
        plotting(maxis,ms)
        save_calib(maxis,ms)


    else:
        print("Wrong choice entered\n")
        calib(tof_spec)

     
#************************************************************************************************************        
#STEP 1
#functions to pass the arguments with the script to run it

def prelim():

    #***********************************
    #pre defined functions needed to run the prelim() function
    def DoTheAnalysis(RunNumber, HDF5path, dset, MolName, MSType,CalFile=None):
        global res
        res = rdl.getRawCleanData(RunNumber, HDF5path, dset,True)
        #plot_TOF(np.mean(res["TOF"], axis=0), msg=MolName)
        PumpMask="FELShutter_markup"
        ProbeMask="ProbeEnergy_markup"
        TOFMS="TOF"
        print(MSType)
        ms = fa.get1dMassSpectrum(res, PumpMask, ProbeMask, TOFMS, MSType)

        fig, ax = plt.subplots()                          #displaying the plots
        ax.plot(np.arange(0, ms.shape[0]), ms)
        ax.set(xlabel='TOF, arb.units', ylabel='Intensity',
               title='1 D TOF spectrum for '+str(MStype)+' data')
        ax.grid()

        np.savetxt(MolName + "_1D_tof_spect_"+str(MStype)+".dat",np.stack((np.arange(0,ms.shape[0]),ms),axis=-1))        #saving the data in text file

        plt.show(block=False)
        print("Press any keyboard key to continue!")
        while True:
            if plt.waitforbuttonpress():
                break
        return ms

    def run(args):
        global RunNumber
        RunNumber = args.RunNumber
        HDF5path = args.HDF5path
        global MolName
        MolName = args.MolName 
        global MStype
        MStype = args.MStype
        LasSet = args.LasSet  
     

        dset = bdsm.SettingsForDataExtraction()
        dset.addMultipleChannelsFromConfigFile('test.inp')
        
        if MStype in {'1','1.','only IR','IR only', 'IR'}:
            MStype = 'pump'
            ms = DoTheAnalysis(RunNumber, HDF5path, dset, MolName,'pump')
            return ms
        elif MStype in {'2','2.','only XUV','XUV only', 'XUV'}:
            MStype = 'probe'
            ms = DoTheAnalysis(RunNumber, HDF5path, dset, MolName,'probe')
            return ms
        elif MStype in {'3','3.','only IR-XUV','IR-XUV only', 'IR-XUV'}:
            MStype = 'pump-probe'
            ms = DoTheAnalysis(RunNumber, HDF5path, dset, MolName,'pump-probe')
            return ms
        elif MStype in {'4','4.','all together','all'}:
            MStype = 'all'
            ms = DoTheAnalysis(RunNumber, HDF5path, dset, MolName,None)
            return ms
        elif MStype in {'5','5.','only background','background'}:
            MStype = 'backgound'
            ms = DoTheAnalysis(RunNumber, HDF5path, dset, MolName,MStype)
            return ms
        else:
            print("wrong choice (Please enter the digits from 1 to 5)")
    
    #***********************************   
    #prelim() function starts here
    parser=argparse.ArgumentParser(description="Step 1\nTakes the information of data to be analyzed !")
    parser.add_argument("-LasSet",help="input laser setting file (default 'test.inp')" ,dest="LasSet", type=str, default="test.inp")
    parser.add_argument("-CalFile",help="input calibration file (default 'CalFile.inp')" ,dest="CalFile", type=str, default="CalFile.inp")
    parser.add_argument("-RunNumber",help="Enter the run number" ,dest="RunNumber", type=int, required=True)
    parser.add_argument("-HDF5path",help="Enter the path of the data" ,dest="HDF5path", type=str, required=True)
    parser.add_argument("-MolName",help="name of the molecule" ,dest="MolName", type=str, required=True)
    parser.add_argument("-Setting",help="For which setting do you want to analyse data?\n1. only IR\n2. only XUV\n3. only IR-XUV\n4. all together\n5. only background.\n Enter the number(1 to 5)" ,dest="MStype", type=str, required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    ms=args.func(args)
    parser=argparse.ArgumentParser()
    return ms


         

#************************************************************************************************************
#Function to show the three main working steps of the script 

def steps(): 

    print("\n-------------------------------------------------------------------------------")
    print("\nBelow are the steps which will be followed for the analysis")
    print("\n-------------------------------------------------------------------------------")

    print("\n1.Specify the following information for the data you want to analyse!")
    print("\n\t- Run Number\n\t- HDF5 File path\n\t- Name of the molecule")
    print("\n\t- Experiment type : IR only/ XUV only/ IR-XUV only/ all/ only background\n")

    print("\n2. Do the calibration of the 1D TOF spectrum ")
    print("(You can add points and save in the file later or directly use the pre created file (CalFile)\n")

    print("\n3. Choose for further analysis\n\t")
    print("- 2D TOF spectrum\n\t- TOF-TOF Covariance Maps\n\t- Ion Yeild v/s pump-probe delay curves\n")


#***********************************************************************************************************
#The program starts here with the main function
# 1. function "steps()" will show the steps involved in the analysis via this script on the console
# 2. function "prelim()" will show the TOF spectrum for the required run number and HDF5 file
# 3. function "calib(ms)" will take the ms in tof values as the x-axis and return the mass calibrated spectrum
# 4. function "MenuAnalysis()" will plot the graphs for 2D mass-spec, pump-probe delay yields, covariance plots 
if __name__=="__main__":
    steps()
    ms= prelim() #Step 1 
    calib(ms) #Step 2
    MenuAnalysis() #Step 3
    










