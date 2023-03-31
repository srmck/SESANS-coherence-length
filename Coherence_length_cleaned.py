from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

def quickpol(rnum):
    w1=Load(str(rnum),LoadMonitors=1)
    MaskDetectors('w1',MaskedWorkspace='MaskWorkspace')
    w1lam=ConvertUnits('w1','Wavelength',AlignBins=1)
    w1m1=ExtractSingleSpectrum('w1_monitors',0)
    w1m1Lam=ConvertUnits('w1m1','Wavelength')
    w1lam=Rebin('w1lam','1.0,0.2,12.0',PreserveEvents=0)
    w1m1Lam=Rebin('w1m1Lam','1.0,0.2,12.0')
    w1norm=w1lam/w1m1Lam
    w1lamInt=SumSpectra('w1norm')
    pol=-1.0*(mtd['w1lamInt_1']-mtd['w1lamInt_2'])/(mtd['w1lamInt_1']+mtd['w1lamInt_2'])
    RenameWorkspace('pol',str(rnum)+'_pol')

def quickpolAlanis(rnum,binning='1.0,0.2,12.0',topandbottom=True):
    w1=Load(str(rnum),LoadMonitors=1)
    MoveInstrumentComponent('w1','SEMSANSWLSFDetector',X=0.5,RelativePosition=True)
    w1=CropWorkspace('w1',StartWorkspaceIndex=40960,EndWorkspaceIndex=40960+63)
    w1lam=ConvertUnits('w1','Wavelength',AlignBins=1)
    w1m1=ExtractSingleSpectrum('w1_monitors',0)
    w1m1Lam=ConvertUnits('w1m1','Wavelength')
    w1lam=Rebin('w1lam',binning,PreserveEvents=0)
    w1m1Lam=Rebin('w1m1Lam',binning)
    w1norm=w1lam/w1m1Lam
    # subtract a per pixel background from the edges of the detector
    wbkgd1=SumSpectra('w1norm',0,10)
    wbkgd=wbkgd1/11.0
    if topandbottom:
        wbkgd2=SumSpectra('w1norm',53,63)
        wbkgd=(wbkgd1+wbkgd2)/22.0
    w1norm=w1norm-wbkgd
    polAll=-1.0*(mtd['w1norm_1']-mtd['w1norm_2'])/(mtd['w1norm_1']+mtd['w1norm_2'])
    RenameWorkspace('polAll',str(rnum)+'_polAll')
    #w1lamInt=SumSpectra('w1norm',20,43)
    # decrease the range of integration to allow for the reduced range of polarised
    # beam on the detector at long wavelengths
    w1lamInt=SumSpectra('w1norm',24,37)
    pol=-1.0*(mtd['w1lamInt_1']-mtd['w1lamInt_2'])/(mtd['w1lamInt_1']+mtd['w1lamInt_2'])
    RenameWorkspace('pol',str(rnum)+'_pol')

def patterson(wstemp, const):
    """Convert workspace into spin echo form"""
    temp = CloneWorkspace(wstemp)
    
    x = temp.extractX()
    x = (x[:, 1:]+x[:, :-1])/2
    temp = Logarithm(temp)
    for i in range(x.shape[0]):
        temp.setY(i, temp.readY(i)/x[i]**2/(temp.sample().getThickness()*0.1))
        temp.setE(i, temp.readE(i)/x[i]**2/(temp.sample().getThickness()*0.1))
        #print('Thickness='+str(temp.sample().getThickness()))
        temp = ConvertUnits(temp, "SpinEchoLength", EFixed=const)
    return temp

def echo_cal2MHz(angle):
    # September 2022 Calibration using GR23
    return 1e3*np.polyval([7.42468285310946e-09, -2.37424465512104e-6, 2.92814364848172e-04,
                       -1.85712645099827e-02, 5.44463590018381e-01], np.abs(angle))
def echo_cal3MHz(angle):
    # September 2022 Calibration using GR23 for 2MHz 
    return 1.5*1e3*np.polyval([7.42468285310946e-09, -2.37424465512104e-6, 2.92814364848172e-04,
                       -1.85712645099827e-02, 5.44463590018381e-01], np.abs(angle))
                       
print('Spin echo constants (-40deg):')
a_2MHz_40deg = echo_cal2MHz(-40)
a_3MHz_40deg = echo_cal3MHz(-40)
print(f'2MHz: {round(a_2MHz_40deg,3)} nm/A^2')    
print(f'3MHz: {round(a_3MHz_40deg,3)} nm/A^2')     

def replotEchoScan(rnum):
    w1=Load(str(rnum),LoadMonitors=1)
    MoveInstrumentComponent('w1','SEMSANSWLSFDetector',X=0.5,RelativePosition=True)
    w1=CropWorkspace('w1',StartWorkspaceIndex=1,EndWorkspaceIndex=1)
    w1lam=ConvertUnits('w1','Wavelength',AlignBins=1)
    w1m1=ExtractSingleSpectrum('w1_monitors',0)
    w1m1Lam=ConvertUnits('w1m1','Wavelength')
    w1lam=Rebin('w1lam','3.0,3.0,12.0',PreserveEvents=0)
    w1m1Lam=Rebin('w1m1Lam','3.0,3.0,12.0')
    w1lamInt=SumSpectra('w1lam')
    w1norm=w1lamInt/w1m1Lam
    wlist=[]
    for  i in range(int(len(mtd['w1norm'].getNames())/2)):
        pol=-1.0*(mtd['w1norm_'+str(i*2+1)]-mtd['w1norm_'+str((i*2)+2)])/(mtd['w1norm_'+str((i*2+1))]+mtd['w1norm_'+str((i*2)+2)])
        RenameWorkspace('pol','pol_'+str(i))
        wlist.append('pol_'+str(i))
    GroupWorkspaces(wlist,OutputWorkspace='pol')
    nper=len(mtd['pol'].getNames())
    xvals=[]
    yvals=[]
    evals=[]
    for j in range(3):
        for i in range(len(mtd['pol'].getNames())):
            xvals.append(i)
            yvals.append(mtd['pol_'+str(i)].dataY(0)[j])
            evals.append(mtd['pol_'+str(i)].dataE(0)[j])
    CreateWorkspace(xvals,yvals,evals,NSpec=3,OutputWorkspace=str(rnum)+'_EchoScan')        

def reduceSESANS2MHz(P0,Sample,PSAng,reload=False,binning='2.5,0.05,12.0'):
    if not mtd.doesExist(str(P0)+'_pol') or reload:
        quickpolAlanis(P0,binning=binning)
    quickpolAlanis(Sample,binning=binning)
    Divide(str(Sample)+'_pol',str(P0)+'_pol',OutputWorkspace=str(Sample)+'_pnorm')
    Divide(str(Sample)+'_polAll',str(P0)+'_polAll',OutputWorkspace=str(Sample)+'_pnormAll')
    temp2=ConvertUnits(str(Sample)+'_pnorm', "SpinEchoLength", EFixed=echo_cal2MHz(PSAng))
    temp3=ConvertUnits(str(Sample)+'_pnormAll', "SpinEchoLength", EFixed=echo_cal2MHz(PSAng))
    temp=patterson(str(Sample)+'_pnorm', echo_cal2MHz(PSAng))
    RenameWorkspace('temp',str(Sample)+'_sesans')
    RenameWorkspace('temp2',str(Sample)+'_PnormSE')
    RenameWorkspace('temp3',str(Sample)+'_PnormAllSE')

def reduceSESANS3MHz(P0,Sample,PSAng,reload=False,binning='5,0.05,12.0',topandbottom=False):
    if not mtd.doesExist(str(P0)+'_pol') or reload:
        quickpolAlanis(P0,binning=binning,topandbottom=topandbottom)
    quickpolAlanis(Sample,binning=binning,topandbottom=topandbottom)
    Divide(str(Sample)+'_pol',str(P0)+'_pol',OutputWorkspace=str(Sample)+'_pnorm')
    Divide(str(Sample)+'_polAll',str(P0)+'_polAll',OutputWorkspace=str(Sample)+'_pnormAll')
    temp2=ConvertUnits(str(Sample)+'_pnorm', "SpinEchoLength", EFixed=echo_cal3MHz(PSAng))
    temp3=ConvertUnits(str(Sample)+'_pnormAll', "SpinEchoLength", EFixed=echo_cal3MHz(PSAng))
    temp=patterson(str(Sample)+'_pnorm', echo_cal3MHz(PSAng))
    RenameWorkspace('temp',str(Sample)+'_sesans')
    RenameWorkspace('temp2',str(Sample)+'_PnormSE')
    RenameWorkspace('temp3',str(Sample)+'_PnormAllSE')
    
def replotEchoScan(rnum):
    w1=Load(str(rnum),LoadMonitors=1)
    MoveInstrumentComponent('w1','SEMSANSWLSFDetector',X=0.5,RelativePosition=True)
    w1=CropWorkspace('w1',StartWorkspaceIndex=1,EndWorkspaceIndex=1)
    w1lam=ConvertUnits('w1','Wavelength',AlignBins=1)
    w1m1=ExtractSingleSpectrum('w1_monitors',0)
    w1m1Lam=ConvertUnits('w1m1','Wavelength')
    w1lam=Rebin('w1lam','3.0,3.0,12.0',PreserveEvents=0)
    w1m1Lam=Rebin('w1m1Lam','3.0,3.0,12.0')
    w1lamInt=SumSpectra('w1lam')
    w1norm=w1lamInt/w1m1Lam
    wlist=[]
    for  i in range(int(len(mtd['w1norm'].getNames())/2)):
        pol=-1.0*(mtd['w1norm_'+str(i*2+1)]-mtd['w1norm_'+str((i*2)+2)])/(mtd['w1norm_'+str((i*2+1))]+mtd['w1norm_'+str((i*2)+2)])
        RenameWorkspace('pol','pol_'+str(i))
        wlist.append('pol_'+str(i))
    GroupWorkspaces(wlist,OutputWorkspace='pol')
    nper=len(mtd['pol'].getNames())
    xvals=[]
    yvals=[]
    evals=[]
    for j in range(3):
        for i in range(len(mtd['pol'].getNames())):
            xvals.append(i)
            yvals.append(mtd['pol_'+str(i)].dataY(0)[j])
            evals.append(mtd['pol_'+str(i)].dataE(0)[j])
    CreateWorkspace(xvals,yvals,evals,NSpec=3,OutputWorkspace=str(rnum)+'_EchoScan')  

def calcAverage(rnums):
    ndiv=len(rnums)
    waverage1=mtd[str(rnums[0])+'_PnormSE']
    waverageAll1=mtd[str(rnums[0])+'_PnormAllSE']
    for i in range(1,len(rnums)):
        waverage1=waverage1+mtd[str(rnums[i])+'_PnormSE']
        waverageAll1=waverageAll1+mtd[str(rnums[i])+'_PnormAllSE']
    waverage=waverage1/float(ndiv)
    waverageAll=waverageAll1/float(ndiv)

def calcAveragePol(rnums):
    ndiv=len(rnums)
    waverage1=mtd[str(rnums[0])+'_pol']
    for i in range(1,len(rnums)):
        waverage1=waverage1+mtd[str(rnums[i])+'_pol']
    waverage_pol=waverage1/float(ndiv)
    
#======================================================================================
#Grating 27: Lines vertical (90 degree orientation) at 2MHz and -40deg
#======================================================================================
binning_2MHz_90deg='3.0,0.05,12.5'
reduceSESANS2MHz(72623,72624,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72625,72626,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72627,72628,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72632,72633,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72634,72635,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72636,72637,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72639,72640,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72641,72642,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72643,72644,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72646,72647,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72648,72649,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72650,72651,-40,True,binning_2MHz_90deg)
reduceSESANS2MHz(72653,72654,-40,True,binning_2MHz_90deg)

sruns=[72624,72626,72628,72633,72635,72637,72640,72642,72644,72647,72649,72651,72654]
calcAverage(sruns)
RenameWorkspace('waverage','pols_2Mhz_90deg')
SaveAscii('pols_2Mhz_90deg',\
'C:\\Users\\xsm\\Documents\\Mantid data\\LARMOR Coherence Feb 2023\\Processed data\\pols_2Mhz_90deg.dat')

#======================================================================================
# Rotate Grating by 82deg so that the lines are at 8deg to the horizontal at 2MHz
#======================================================================================
# quick check run of 2uamps each
binning1='3.0,0.05,12.5'

reduceSESANS2MHz(72665,72666,-40,True,binning1)
reduceSESANS2MHz(72667,72668,-40,True,binning1)
reduceSESANS2MHz(72669,72670,-40,True,binning1)
reduceSESANS2MHz(72672,72673,-40,True,binning1)
reduceSESANS2MHz(72674,72675,-40,True,binning1)
reduceSESANS2MHz(72674,72675,-40,True,binning1)
reduceSESANS2MHz(72676,72677,-40,True,binning1)
reduceSESANS2MHz(72679,72680,-40,True,binning1)
reduceSESANS2MHz(72681,72682,-40,True,binning1)
reduceSESANS2MHz(72683,72684,-40,True,binning1)
reduceSESANS2MHz(72686,72687,-40,True,binning1)
reduceSESANS2MHz(72688,72689,-40,True,binning1)
reduceSESANS2MHz(72690,72691,-40,True,binning1)
reduceSESANS2MHz(72693,72694,-40,True,binning1)
reduceSESANS2MHz(72695,72696,-40,True,binning1)
reduceSESANS2MHz(72697,72698,-40,True,binning1)
reduceSESANS2MHz(72700,72701,-40,True,binning1)
reduceSESANS2MHz(72702,72703,-40,True,binning1)



sruns=[72666,72668,72670,72673,72675,72677,72680,72682,72684,72687,72689,72691,72694,72696,72698,72701,72703]
calcAverage(sruns)
RenameWorkspace('waverage','waverage8deg')

base_2mhz=mtd['waverage8deg']*1.0
xvals=base_2mhz.dataX(0)
for i in range(len(xvals)-1):
    x1=0.5*(xvals[i]+xvals[i+1])
    base_2mhz.dataY(0)[i]=0.33+0.67*np.cos(0.57*np.sqrt(x1/1000.0))
    base_2mhz.dataE(0)[i]=0.0

base_2mhz_alter=mtd['waverage8deg']*1.0
xvals=base_2mhz_alter.dataX(0)
for i in range(len(xvals)-1):
    x1=0.5*(xvals[i]+xvals[i+1])
    base_2mhz_alter.dataY(0)[i]=0.35+0.65*np.cos(0.57*np.sqrt(x1/1000.0))
    base_2mhz_alter.dataE(0)[i]=0.0





#waverage8deg_bkgd=mtd['waverage8deg']-wbase8deg


# New function after Sam's fit
wbase8deg_new=mtd['waverage8deg']*1.0

xvals=wbase8deg_new.dataX(0)
for i in range(len(xvals)-1):
    x1=0.5*(xvals[i]+xvals[i+1])
    wbase8deg_new.dataY(0)[i]=0.33+0.67*np.cos(0.57*np.sqrt(x1/1000.0))
    wbase8deg_new.dataE(0)[i]=0.0
waverage8deg_bkgd_new=mtd['waverage8deg']-wbase8deg_new

#======================================================================================
# Grating still at 8deg to the horizontal and moved to 3MHz
#======================================================================================
#Quick runs to check stability
binning1='5.5,0.1,13.3'

reduceSESANS3MHz(72707,72708,-40,True,binning1)
reduceSESANS3MHz(72709,72708,-40,True,binning1)

quickpolAlanis(72712,binning1)
quickpolAlanis(72713,binning1)
quickpolAlanis(72714,binning1)

# After waiting for 3MHz to stabilise
reduceSESANS3MHz(72716,72717,-40,True,binning1)
reduceSESANS3MHz(72718,72719,-40,True,binning1)
reduceSESANS3MHz(72720,72721,-40,True,binning1)

reduceSESANS3MHz(72723,72724,-40,True,binning1)
reduceSESANS3MHz(72725,72726,-40,True,binning1)
reduceSESANS3MHz(72727,72728,-40,True,binning1)
reduceSESANS3MHz(72730,72731,-40,True,binning1)
reduceSESANS3MHz(72732,72733,-40,True,binning1)
reduceSESANS3MHz(72734,72735,-40,True,binning1)
reduceSESANS3MHz(72737,72738,-40,True,binning1)
reduceSESANS3MHz(72739,72740,-40,True,binning1)
reduceSESANS3MHz(72741,72742,-40,True,binning1)
reduceSESANS3MHz(72744,72745,-40,True,binning1)
reduceSESANS3MHz(72746,72747,-40,True,binning1)
reduceSESANS3MHz(72748,72749,-40,True,binning1)
reduceSESANS3MHz(72751,72752,-40,True,binning1)

sruns=[72717,72719,72721,72724,72726,72728,72731,72733,72735,72738,72740,72742,72745,72747,72749,72752]
calcAverage(sruns)
RenameWorkspace('waverage','waverage8deg3mhz')

#Create cos theta background for 3MHz data 
base_3mhz=mtd['waverage8deg3mhz']*1.0

xvals=base_3mhz.dataX(0)
for i in range(len(xvals)-1):
    x1=0.5*(xvals[i]+xvals[i+1])
    base_3mhz.dataY(0)[i]=0.35+0.65*np.cos(0.57*np.sqrt(x1/1000.0*2/3))
    base_3mhz.dataE(0)[i]=0.0

base_3mhz_alter=mtd['waverage8deg3mhz']*1.0
xvals=base_3mhz_alter.dataX(0)
for i in range(len(xvals)-1):
    x1=0.5*(xvals[i]+xvals[i+1])
    base_3mhz_alter.dataY(0)[i]=0.33+0.67*np.cos(0.57*np.sqrt(x1/1000.0*2/3))
    base_3mhz_alter.dataE(0)[i]=0.0
Gr27_3mhz_8deg_bgnd_sub=mtd['waverage8deg3mhz']-base_3mhz_alter


#waverage8deg3mhz_bkgd_sub=mtd['waverage8deg3mhz']-waverage8deg3mhz_bgnd


#Write data to ascii
SaveAscii('waverage8deg3mhz',r'U:\Users\Pynn\Feb_2023\GR27_8deg_3Mhz0p5_bin_bgnd_sub.txt')
#SaveRKH('GR27_Average_0p05_bin',r'U:\Users\Pynn\Feb_2023\GR27_Average_0p5_bin.txt',Append=False)

#======================================================================================
# Grating 27 at 5deg to the horizontal and moved to 3MHz
#======================================================================================
binning1='5.5,0.1,13.3'

reduceSESANS3MHz(72754,72755,-40,True,binning1)
sruns=[72755]
calcAverage(sruns)
RenameWorkspace('waverage','waverage5deg3mhz')
#======================================================================================
# Grating 27 at 4p5deg to the horizontal and moved to 3MHz
#======================================================================================
binning1='5.5,0.1,13.3'

#quickpol(72795)

reduceSESANS3MHz(72757,72758,-40,True,binning1) # Quick run   12 microamps

reduceSESANS3MHz(72761,72762,-40,True,binning1)  # Long run 20 microamps
reduceSESANS3MHz(72763,72764,-40,True,binning1)
reduceSESANS3MHz(72765,72766,-40,True,binning1)
reduceSESANS3MHz(72768,72769,-40,True,binning1)
reduceSESANS3MHz(72770,72771,-40,True,binning1)
reduceSESANS3MHz(72772,72773,-40,True,binning1)
reduceSESANS3MHz(72775,72776,-40,True,binning1)
reduceSESANS3MHz(72777,72778,-40,True,binning1)
reduceSESANS3MHz(72779,72780,-40,True,binning1)
reduceSESANS3MHz(72782,72783,-40,True,binning1)
reduceSESANS3MHz(72784,72785,-40,True,binning1)
reduceSESANS3MHz(72786,72787,-40,True,binning1)
reduceSESANS3MHz(72789,72790,-40,True,binning1)
reduceSESANS3MHz(72791,72792,-40,True,binning1)
reduceSESANS3MHz(72793,72794,-40,True,binning1)
reduceSESANS3MHz(72796,72797,-40,True,binning1)
reduceSESANS3MHz(72798,72799,-40,True,binning1)
reduceSESANS3MHz(72800,72801,-40,True,binning1)
reduceSESANS3MHz(72803,72804,-40,True,binning1)
reduceSESANS3MHz(72805,72806,-40,True,binning1)
reduceSESANS3MHz(72807,72808,-40,True,binning1)
reduceSESANS3MHz(72810,72811,-40,True,binning1)
reduceSESANS3MHz(72812,72813,-40,True,binning1)
reduceSESANS3MHz(72814,72815,-40,True,binning1)
reduceSESANS3MHz(72817,72818,-40,True,binning1)

sruns=[72762,72764,72766,72769,72771,72773,72776,72778,72780,72783,72785,72787,72790,72792,72794,72797,72799,72801,72804,72806,72808,72811,72813,72815,72818]
calcAverage(sruns)
RenameWorkspace('waverage','waverage4p5deg3mhz')

SaveAscii('waverage4p5deg3mhz',r'U:\Users\Pynn\Feb_2023\waverage4p5deg3mhz.txt')

#======================================================================================
# Grating 27 at 90deg to the horizontal at 3MHz
#======================================================================================
binning1='5.5,0.05,13.3'

#reduceSESANS3MHz(72820,72821,-40,True,binning1) # Quick run 2 microamps

reduceSESANS3MHz(72823,72824,-40,True,binning1)  # Long run 20 microamps
reduceSESANS3MHz(72825,72826,-40,True,binning1)
reduceSESANS3MHz(72827,72828,-40,True,binning1)
reduceSESANS3MHz(72830,72831,-40,True,binning1)
reduceSESANS3MHz(72832,72833,-40,True,binning1)
reduceSESANS3MHz(72834,72835,-40,True,binning1)
reduceSESANS3MHz(72837,72838,-40,True,binning1)
reduceSESANS3MHz(72839,72840,-40,True,binning1)
reduceSESANS3MHz(72841,72842,-40,True,binning1)
reduceSESANS3MHz(72844,72845,-40,True,binning1)
reduceSESANS3MHz(72846,72847,-40,True,binning1)
reduceSESANS3MHz(72848,72849,-40,True,binning1)
reduceSESANS3MHz(72851,72852,-40,True,binning1)
reduceSESANS3MHz(72853,72854,-40,True,binning1)
reduceSESANS3MHz(72855,72856,-40,True,binning1)
reduceSESANS3MHz(72858,72859,-40,True,binning1)
reduceSESANS3MHz(72860,72861,-40,True,binning1)
reduceSESANS3MHz(72862,72863,-40,True,binning1)
reduceSESANS3MHz(72865,72866,-40,True,binning1)
reduceSESANS3MHz(72867,72868,-40,True,binning1)
reduceSESANS3MHz(72869,72870,-40,True,binning1)
reduceSESANS3MHz(72872,72873,-40,True,binning1)
reduceSESANS3MHz(72874,72875,-40,True,binning1)
reduceSESANS3MHz(72876,72877,-40,True,binning1)
reduceSESANS3MHz(72879,72880,-40,True,binning1)
reduceSESANS3MHz(72881,72882,-40,True,binning1)
reduceSESANS3MHz(72883,72884,-40,True,binning1)
reduceSESANS3MHz(72886,72887,-40,True,binning1)

sruns=[72824,72826,72828,72831,72833,72835,72838,72840,72842,72845,72847,72849,72852,72854,72856,72859,72861,\
72863,72866,72868,72870,72873,72875,72877,72880,72882,72884,72887]
calcAverage(sruns)
RenameWorkspace('waverage','waverage90deg3mhz')

SaveAscii('waverage90deg3mhz',r'U:\Users\Pynn\Feb_2023\waverage90deg3mhz.txt')

binning1='5.5,0.025,13.3'

#reduceSESANS3MHz(72820,72821,-40,True,binning1) # Quick run 2 microamps

reduceSESANS3MHz(72823,72824,-40,True,binning1)  # Long run 20 microamps
reduceSESANS3MHz(72825,72826,-40,True,binning1)
reduceSESANS3MHz(72827,72828,-40,True,binning1)
reduceSESANS3MHz(72830,72831,-40,True,binning1)
reduceSESANS3MHz(72832,72833,-40,True,binning1)
reduceSESANS3MHz(72834,72835,-40,True,binning1)
reduceSESANS3MHz(72837,72838,-40,True,binning1)
reduceSESANS3MHz(72839,72840,-40,True,binning1)
reduceSESANS3MHz(72841,72842,-40,True,binning1)
reduceSESANS3MHz(72844,72845,-40,True,binning1)
reduceSESANS3MHz(72846,72847,-40,True,binning1)
reduceSESANS3MHz(72848,72849,-40,True,binning1)
reduceSESANS3MHz(72851,72852,-40,True,binning1)
reduceSESANS3MHz(72853,72854,-40,True,binning1)
reduceSESANS3MHz(72855,72856,-40,True,binning1)
reduceSESANS3MHz(72858,72859,-40,True,binning1)
reduceSESANS3MHz(72860,72861,-40,True,binning1)
reduceSESANS3MHz(72862,72863,-40,True,binning1)
reduceSESANS3MHz(72865,72866,-40,True,binning1)
reduceSESANS3MHz(72867,72868,-40,True,binning1)
reduceSESANS3MHz(72869,72870,-40,True,binning1)
reduceSESANS3MHz(72872,72873,-40,True,binning1)
reduceSESANS3MHz(72874,72875,-40,True,binning1)
reduceSESANS3MHz(72876,72877,-40,True,binning1)
reduceSESANS3MHz(72879,72880,-40,True,binning1)
reduceSESANS3MHz(72881,72882,-40,True,binning1)
reduceSESANS3MHz(72883,72884,-40,True,binning1)
reduceSESANS3MHz(72886,72887,-40,True,binning1)

sruns=[72824,72826,72828,72831,72833,72835,72838,72840,72842,72845,72847,72849,72852,72854,72856,72859,72861,\
72863,72866,72868,72870,72873,72875,72877,72880,72882,72884,72887]
calcAverage(sruns)
RenameWorkspace('waverage','waverage90deg3mhz_0p025')


