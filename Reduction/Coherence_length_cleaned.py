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
binning_2MHz_90deg='3.0,0.05,12.5' #original binning

empties = [72623,72625,72627,72632,72634,72636,72639,72641,72643,72646,72648,72650,72653,72655]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS2MHz(r,r+1,-40,True,binning=binning_2MHz_90deg)
calcAverage(samples)
data_name = 'pols_2Mhz_90deg'+str(binning_2MHz_90deg).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_2MHz_90deg_2='3.0,0.05,13.5'

empties = [72623,72625,72627,72632,72634,72636,72639,72641,72643,72646,72648,72650,72653,72655]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS2MHz(r,r+1,-40,True,binning=binning_2MHz_90deg_2)
calcAverage(samples)
data_name = 'pols_2Mhz_90deg'+str(binning_2MHz_90deg_2).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_2MHz_90deg_3='3.0,0.025,13.5'

empties = [72623,72625,72627,72632,72634,72636,72639,72641,72643,72646,72648,72650,72653,72655]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS2MHz(r,r+1,-40,True,binning=binning_2MHz_90deg_3)
calcAverage(samples)
data_name = 'pols_2Mhz_90deg'+str(binning_2MHz_90deg_3).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')
#======================================================================================
# Rotate Grating by 82deg so that the lines are at 8deg to the horizontal at 2MHz
#======================================================================================
binning_2MHz_8deg='3.0,0.05,12.5'

empties = [72665,72667,72669,72672,72674,72676,72679,72681,72683,72686,72688,72690,72693,72695,72697,72700,72702]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS2MHz(r,r+1,-40,True,binning=binning_2MHz_8deg)
calcAverage(samples)
data_name = 'pols_2Mhz_8deg'+str(binning_2MHz_8deg).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_2MHz_8deg_2='3.0,0.025,13.5'

empties = [72665,72667,72669,72672,72674,72676,72679,72681,72683,72686,72688,72690,72693,72695,72697,72700,72702]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS2MHz(r,r+1,-40,True,binning=binning_2MHz_8deg_2)
calcAverage(samples)
data_name = 'pols_2Mhz_8deg'+str(binning_2MHz_8deg_2).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_2MHz_8deg_3='3.0,0.05,13.5'

empties = [72665,72667,72669,72672,72674,72676,72679,72681,72683,72686,72688,72690,72693,72695,72697,72700,72702]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS2MHz(r,r+1,-40,True,binning=binning_2MHz_8deg_3)
calcAverage(samples)
data_name = 'pols_2Mhz_8deg'+str(binning_2MHz_8deg_3).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

#======================================================================================
# Grating still at 8deg to the horizontal; rf frequency increased to 3MHz
#======================================================================================
binning_3MHz_8deg='5.5,0.05,13.25'

empties = [72716,72718,72720,72723,72725,72727,72730,72732,72734,72737,72739,72741,72744,72746,72748,72751]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_8deg)
calcAverage(samples)
data_name = 'pols_3Mhz_8deg'+str(binning_3MHz_8deg).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_3MHz_8deg_2='5.5,0.075,13.25'

empties = [72716,72718,72720,72723,72725,72727,72730,72732,72734,72737,72739,72741,72744,72746,72748,72751]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_8deg_2)
calcAverage(samples)
data_name = 'pols_3Mhz_8deg'+str(binning_3MHz_8deg_2).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_3MHz_8deg_3='5.5,0.1,13.25'

empties = [72716,72718,72720,72723,72725,72727,72730,72732,72734,72737,72739,72741,72744,72746,72748,72751]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_8deg_3)
calcAverage(samples)
data_name = 'pols_3Mhz_8deg'+str(binning_3MHz_8deg_3).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

#======================================================================================
# Grating 27 at 4p5deg to the horizontal; rf frequency at 3MHz
#======================================================================================
binning_3MHz_4p5deg='5.5,0.1,13.25'

empties = [72761,72763,72765,72768,72770,72772,72775,72777,72779,72782,72784,72786,72789,72791,72793,72796,72798,72800,72803,72805,72807,72810,72812,72814,72817]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_4p5deg)
calcAverage(samples)
data_name = 'pols_3Mhz_4p5deg'+str(binning_3MHz_4p5deg).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_3MHz_4p5deg_2='5.5,0.05,13.25'

empties = [72761,72763,72765,72768,72770,72772,72775,72777,72779,72782,72784,72786,72789,72791,72793,72796,72798,72800,72803,72805,72807,72810,72812,72814,72817]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_4p5deg_2)
calcAverage(samples)
data_name = 'pols_3Mhz_4p5deg'+str(binning_3MHz_4p5deg_2).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

#======================================================================================
# Grating 27 at 90deg to the horizontal; rf frequency at 3MHz
#======================================================================================
binning_3MHz_90deg='5.5,0.05,13.25'

empties = [72823,72825,72827,72830,72832,72834,72837,72839,72841,72844,72846,72848,72851,72853,72855,72858,72860,72862,72865,72867,72869,72872,72874,72876,72879,72881,72883,72886,72888]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_90deg)
calcAverage(samples)
data_name = 'pols_3Mhz_90deg'+str(binning_3MHz_90deg).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_3MHz_90deg_1='5.5,0.1,13.25'

empties = [72823,72825,72827,72830,72832,72834,72837,72839,72841,72844,72846,72848,72851,72853,72855,72858,72860,72862,72865,72867,72869,72872,72874,72876,72879,72881,72883,72886,72888]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_90deg_1)
calcAverage(samples)
data_name = 'pols_3Mhz_90deg'+str(binning_3MHz_90deg_1).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')

binning_3MHz_90deg_3='5.5,0.025,13.25'

empties = [72823,72825,72827,72830,72832,72834,72837,72839,72841,72844,72846,72848,72851,72853,72855,72858,72860,72862,72865,72867,72869,72872,72874,72876,72879,72881,72883,72886,72888]
samples = [e+1 for e in empties]
for r in empties:
    reduceSESANS3MHz(r,r+1,-40,True,binning=binning_3MHz_90deg_3)
calcAverage(samples)
data_name = 'pols_3Mhz_90deg'+str(binning_3MHz_90deg_3).replace(',','_')
RenameWorkspace('waverage',data_name)
SaveAscii(data_name,\
'C:\\Users\\xsm\\Documents\\\GitHub\\SESANS-coherence-length\\Processed data\\'+data_name+'.dat')




