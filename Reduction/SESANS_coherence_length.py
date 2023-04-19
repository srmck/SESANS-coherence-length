import sys
from genie_python import genie as g, BLOCK_NAMES as b
sys.path.append(r"c:\instrument\scripts")
from instrument.larmor import * # pylint: disable=wildcard-import, unused-wildcard-import
import LSS.SESANSroutines as ss
import time as time

def setdaeAlanisScanning():
    g.change_start()
    g.change_tables(detector=r"C:\Instrument\Settings\config\NDXLARMOR\configurations\tables\Alanis_Detector.dat")
    g.change_tables(spectra=r"C:\Instrument\Settings\config\NDXLARMOR\configurations\tables\spectra_scanning_Alanis.dat")
    g.change_tables(wiring=r"C:\Instrument\Settings\config\NDXLARMOR\configurations\tables\Alanis_Wiring_dae3.dat")
    g.change_tcb(low=5.0,high=100000.0,step=100.0,trange=1,log=0)
    g.change_finish()

def deGaussDCMagnets():
    fields=[15,-10,5,-2.5,1.25,-0.6,0]
    for i in fields:
        g.cset(DCMagField1=i)
        time.sleep(0.5)
        g.cset(DCMagField2=i)
        time.sleep(0.5)
        g.cset(DCMagField3=i)
        time.sleep(0.5)
        g.cset(DCMagField4=i)
        time.sleep(30)


def zerofields():
    g.cset(DCMagField1=0.0)
    time.sleep(0.5)
    g.cset(DCMagField2=0.0)
    time.sleep(0.5)
    g.cset(DCMagField3=0.0)
    time.sleep(0.5)
    g.cset(DCMagField4=0.0)
    

def rezerofields3MHz():
    g.cset(DCMagField1=0.0)
    time.sleep(0.5)
    g.cset(DCMagField2=0.0)
    time.sleep(0.5)
    g.cset(DCMagField3=0.0)
    time.sleep(0.5)
    g.cset(DCMagField4=0.0)
    time.sleep(30)
    g.cset(DCMagField1=103.1)
    time.sleep(0.5)
    g.cset(DCMagField2=103.1)
    time.sleep(0.5)
    g.cset(DCMagField3=-104.0)
    time.sleep(0.5)
    g.cset(DCMagField4=-104.0)
    time.sleep(30)
 
def rezerofields2MHz():
    g.cset(DCMagField1=0.0)
    time.sleep(0.5)
    g.cset(DCMagField2=0.0)
    time.sleep(0.5)
    g.cset(DCMagField3=0.0)
    time.sleep(0.5)
    g.cset(DCMagField4=0.0)
    time.sleep(30)
    g.cset(DCMagField1=69.0)
    time.sleep(0.5)
    g.cset(DCMagField2=69.0)
    time.sleep(0.5)
    g.cset(DCMagField3=-69.6)
    time.sleep(0.5)
    g.cset(DCMagField4=-69.6)
    time.sleep(30)
    
def rezerofields1MHz():
    g.cset(DCMagField1=0.0)
    time.sleep(0.5)
    g.cset(DCMagField2=0.0)
    time.sleep(0.5)
    g.cset(DCMagField3=0.0)
    time.sleep(0.5)
    g.cset(DCMagField4=0.0)
    time.sleep(30)
    g.cset(DCMagField1=34.5)
    time.sleep(0.5)
    g.cset(DCMagField2=34.5)
    time.sleep(0.5)
    g.cset(DCMagField3=-34.8)
    time.sleep(0.5)
    g.cset(DCMagField4=-34.8)
    time.sleep(30) 

def rezerofields0p5MHz():
    g.cset(DCMagField1=0.0)
    time.sleep(0.5)
    g.cset(DCMagField2=0.0)
    time.sleep(0.5)
    g.cset(DCMagField3=0.0)
    time.sleep(0.5)
    g.cset(DCMagField4=0.0)
    time.sleep(30)
    g.cset(DCMagField1=17.25)
    time.sleep(0.5)
    g.cset(DCMagField2=17.25)
    time.sleep(0.5)
    g.cset(DCMagField3=-17.4)
    time.sleep(0.5)
    g.cset(DCMagField4=-17.4)
    time.sleep(30)    


def TuesMorningLong():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        setup_dae_event()
        setup_dae_alanis() 
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        #ss.set_poleshoe_angle2(-90,1177.2,2)
        #rezerofields2MHz()    
        
        #-40deg
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        #ss.set_poleshoe_angle2(-40,1177.2,2)
        ss.auto_tune("Echo_Coil_SP", 6300, 7300, 21, 50, "Echo Scan 2MHz -40deg", True)
        setup_dae_event()
        setup_dae_alanis()        
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')


def TuesAfternoonLong():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        #ss.set_poleshoe_angle2(-90,1177.2,2)
        #rezerofields2MHz()    
        
        #-40deg
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        #ss.set_poleshoe_angle2(-40,1177.2,2)
        ss.auto_tune("Echo_Coil_SP", 6300, 7300, 21, 50, "Echo Scan 2MHz -40deg", True)
        setup_dae_event()
        setup_dae_alanis()        
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 90deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
def WedMorning():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        #ss.set_poleshoe_angle2(-90,1177.2,2)
        #rezerofields2MHz()    
        
        #-40deg
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        #ss.set_poleshoe_angle2(-40,1177.2,2)
        ss.auto_tune("Echo_Coil_SP", 6300, 7300, 21, 50, "Echo Scan 2MHz -40deg", True)
        setup_dae_event()
        setup_dae_alanis()        
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 8deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 8deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 8deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 8deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        do_sans(title='{Empty Beam a1=14x14 a2=4x4 8deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

        g.cset(changertranslation=236) #center of grating
        g.waitfor_move()
        do_sans(title='{GR27 8deg 2MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

def ThursdayMorning():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        ss.auto_tune("Echo_Coil_SP", 7000, 7700, 21, 50, "Echo Scan 3MHz -40deg", True)

        setup_dae_event()
        setup_dae_alanis()        
        
        for i in range(3):
            g.cset(changertranslation=236-70) #empty beam
            g.waitfor_move()
            do_sans(title='{Empty Beam a1=14x14 a2=4x4 3MHz -40poleshoe}', uamps=10, thickness=1, dae='sesans')

            #g.cset(changertranslation=236) #center of grating
            #g.waitfor_move()
            #do_sans(title='{GR27 8deg 3MHz -40poleshoe}', uamps=10, thickness=1, dae='sesans')
            
def ThursdayAfternoon():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        ss.auto_tune("Echo_Coil_SP", 7000, 7700, 21, 50, "Echo Scan 3MHz -40deg", True)

        setup_dae_event()
        setup_dae_alanis()        
        
        for i in range(3):
            g.cset(changertranslation=236-70) #empty beam
            g.waitfor_move()
            do_sans(title='{Empty Beam a1=14x14 a2=4x4 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

            g.cset(changertranslation=236) #center of grating
            g.waitfor_move()
            do_sans(title='{GR27 8deg 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')            

def FriMorning():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        ss.auto_tune("Echo_Coil_SP", 7000, 7700, 21, 50, "Echo Scan 3MHz -40deg", True)

        setup_dae_event()
        setup_dae_alanis()        
        
        for i in range(3):
            g.cset(changertranslation=236-70) #empty beam
            g.waitfor_move()
            do_sans(title='{Empty Beam a1=14x14 a2=4x4 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

            g.cset(changertranslation=236) #center of grating
            g.waitfor_move()
            do_sans(title='{GR27 5deg 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')         


def FriLunch():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        ss.auto_tune("Echo_Coil_SP", 7000, 7700, 21, 50, "Echo Scan 3MHz -40deg", True)

        setup_dae_event()
        setup_dae_alanis()        
        
        for i in range(3):
            g.cset(changertranslation=236-70) #empty beam
            g.waitfor_move()
            do_sans(title='{Empty Beam a1=14x14 a2=4x4 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

            g.cset(changertranslation=236) #center of grating
            g.waitfor_move()
            do_sans(title='{GR27 4p5deg 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')         

def SatEvening():

    g.change_sample_par('Width', '4')
    g.change_sample_par('height', '4')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:
        
        g.cset(changertranslation=236-70) #empty beam
        g.waitfor_move()
        ss.auto_tune("Echo_Coil_SP", 7100, 7800, 21, 50, "Echo Scan 3MHz -40deg", True)

        setup_dae_event()
        setup_dae_alanis()        
        
        for i in range(3):
            g.cset(changertranslation=236-70) #empty beam
            g.waitfor_move()
            do_sans(title='{Empty Beam a1=14x14 a2=4x4 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')

            g.cset(changertranslation=236) #center of grating
            g.waitfor_move()
            do_sans(title='{GR27 90deg 3MHz -40poleshoe}', uamps=20, thickness=1, dae='sesans')  
            
            
            
def MonMorning_refraction():

    g.change_sample_par('Width', '5')
    g.change_sample_par('height', '5')
    g.change_sample_par('Geometry', 'Flat Plate')
        
    while True:

        setup_dae_event()
        setup_dae_alanis()        
        
        for i in range(3):
            do_trans(title='{Empty Beam a1=14x14 a2=5x5 3MHz -40poleshoe}', pos='1R',uamps=5, thickness=1)
            do_sans(title='{Empty Beam a1=14x14 a2=5x5 3MHz -40poleshoe}', pos='1R',uamps=20, thickness=1, dae='sesans')

            do_trans(title='{Sample2_Al2O3_H2o_4p5m 3MHz -40poleshoe}', pos='2R',uamps=5, thickness=1)            
            do_sans(title='{Sample2_Al2O3_H2o_4p5m 3MHz -40poleshoe}', pos='2R',uamps=20, thickness=1, dae='sesans')            

            do_trans(title='{Sample3_Al2O3_d2o_4p5m 3MHz -40poleshoe}', pos='3R',uamps=5, thickness=1)
            do_sans(title='{Sample3_Al2O3_d2o_4p5m 3MHz -40poleshoe}', pos='3R',uamps=20, thickness=1, dae='sesans')            
            
        g.cset(samplepos='1R')
        g.waitfor_move()
        ss.auto_tune("Echo_Coil_SP", 7100, 7800, 21, 50, "Echo Scan 3MHz -40deg", True)
