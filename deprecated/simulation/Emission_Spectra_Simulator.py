"""

Fifth Attempt at making the emission spectra simulator


"""

import os
import sys
import itertools
import numpy as np
import pandas as pd

from scipy import stats
import scipy.integrate
import scipy.optimize
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Main.atmosphere_effects.biosig_molecule import *

import SEAS_Utils as utils
import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.data_loader import *
from SEAS_Utils.common_utils.constants import *

import SEAS_Aux.cross_section.cross_section_calculator as csc
from SEAS_Utils.common_utils.DIRs import Temp_DIR, HITRAN_Water_Lines,HITRAN_Lines,HITRAN_CH4_Lines



def read_modtran_input():
    
    with open("modtran_atmosphere_profile.txt","r") as f:
        data = f.read().split("\n")
        
        ZArray,PArray,TArray = [],[],[]
        for row in data:
            if "#" in row:
                continue
            Z,P,T = row.split()[1:4]
            
            if float(P) == 0.0:
                P = 0.001
            ZArray.append(float(Z)*1000)
            PArray.append(float(P)*100)
            TArray.append(float(T))
    
    return np.array(ZArray),np.array(PArray),np.array(TArray)

def load_Earth_MR():

    MR_file = "../../input/atmosphere_data/Mixing_Ratio/earth.txt"
    
    data = multi_column_file_loader(MR_file,type="mixed")
    
    molecules, molecule_MR = [],[]
    for i in data[1:]:
        molecules.append(i[0])
        molecule_MR.append([float(x)/100. for x in i[1:]])

    MR_pressure = [float(x) for x in data[0][1:]]
    
    return molecules,molecule_MR,MR_pressure

def interpolate_MR(molecule_MR,MR_pressure,normalized_pressure):
    
    normalized_molecule_MR = []
    x = np.log(MR_pressure)
    X = np.log(normalized_pressure) 
    
    for y in molecule_MR:
        
        interp = interp1d(x,y,fill_value="extrapolate")
        interpolated_MR = interp(X)
        
        normalized_molecule_MR.append(interpolated_MR)
    
    #normalized_abundance = np.array(normalized_molecule_MR).T
    return normalized_molecule_MR  

def load_base_MR(stellar_type = "Sun", 
                 main_gas = "H2",
                 mu_atm = 4.6*amu2g, 
                 g_atm = 3200.):
     
    
    """
    CH4 = "CH4",
    CH4_base_mr = 1e-6,
    CH4_mr_scale_factor = 10.,
    CH4_cutoff = 10,   
    
    generate mixing ratio profile from renyu's simulation
    
    stellar_type = which type of star to consider: Sun, Ma, Mq
    main_gas = dominant composition of the atmosphere
    base_mr = base mixing ratio of biosignature gas in consideration
    mu_atm = average molar mass of atmosphere for each scenario (grams)
    g_atm = gravitational constant of atmosphere for each scenario (cm s**-2)
    """
    datapath = "../../input/concentration"
    tp_filename = "../../input/atmosphere_data/TP_Profile_Calculated/%s_%s.txt"%(main_gas,stellar_type)
    #mr_filename_ref = "../../input/atmosphere_data/Mixing_Ratio_Calculated/%s_%s_ref.txt"%(main_gas,stellar_type)    
    mr_filename_ref = "../../input/atmosphere_data/Mixing_Ratio_Calculated/%s_%s_ref_vari.txt"%(main_gas,stellar_type)
    
    ########################
    ###Read in key file, to translate between Hu code numbering convention and English names.
    ########################
    species_name_file= os.path.join(datapath,'SpeciesName.dat') #Directory storing the key mapping species to the integers 1-111 used to ID them in the code
    imported_key = np.genfromtxt(species_name_file, skip_header=1, skip_footer=0, delimiter='\t', dtype=None) #Import mapping between numerical ID in code and species name.
    
    key_species2num={} #initalize empty dictionary
    key_num2species={} #initalize empty dictionary, to hold mapping FROM species name TO species number
    
    for i in range(0, len(imported_key)):
        species_number, species_name = imported_key[i]
        key_species2num[species_name]=species_number #Holds mapping FROM species name TO species number
        key_num2species[species_number]=species_name #Holds mapping FROM species number TO species name
    
    scenario_file = os.path.join(datapath,'ConcentrationSTD_%s_%s.dat'%(main_gas, stellar_type))
    output_data=np.genfromtxt(scenario_file, skip_header=2, skip_footer=0, unpack=False) 

    z_centers=output_data[:,0]*km2cm # Center of altitude bins, km converted to cm
    z_lower=output_data[:,1]*km2cm # Lower edge of altitude bins, km converted to cm
    z_upper=output_data[:,2]*km2cm # Upper edge of altitude bins, km converted to cm
    T_z=output_data[:,3] # Temperature(z), in K
    P_z=output_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    P_z_Pa = output_data[:,4] # Pressure(z), in Pa 
    n_z_species=output_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    # Saving the TP profile
    newP,newT,z_out = [],[],[]
    for P,T,Z in zip(P_z_Pa,T_z,z_centers/100): # zcenter in meter
        if P < 0.01:
            break
        newT.append(T)
        newP.append(P)
        z_out.append(Z)
    PT_Profile = [newP,newT,z_out]
    
    MR_Profile = {}
    #sum number densities across species. This is a profile for the whole atmosphere.
    n_z=np.sum(n_z_species,1)     

    for j,info in enumerate(n_z_species.T):
        l = key_num2species[j+1]
        if max(info/n_z) > 10**-7:
            MR_Profile[l] = (info/n_z)[:len(newP)]
    
    """
    # add methane into the fold
    z_widths=z_upper-z_lower 
    colden_atm=np.sum(n_z*z_widths)
    colden_isoprene=CH4_base_mr*colden_atm
    k=1.380658e-16 #boltzmann constant, erg/K
    H=k*T_z[0]/(mu_atm*g_atm) #atmospheric scale height
    zs = H*CH4_mr_scale_factor
    k=1.380658e-23 #boltzmann constant, erg/K
    mr_z_unscaled=np.exp(-z_centers/zs)
    norm_const=colden_isoprene/np.sum(mr_z_unscaled*n_z*z_widths)
    mr_z=norm_const*mr_z_unscaled
 
    mr_z = mr_z[:len(newP)]
    MR_Profile[CH4] = mr_z
    """
    
    return PT_Profile, MR_Profile





def load_Earth_PT():

    P,T = np.genfromtxt("../../input/atmosphere_data/TP_profile/earth.txt").T
    
    return P,T

def interpolate_PT_profile(PT_Pressure,PT_Temperature,normalized_pressure):
    
    x = np.log(TP_Pressure)
    y = TP_Temperature
    
    
    x = np.concatenate([[100],x,[-100]])
    y = np.concatenate([[y[0]],y,[y[-1]]])
    
    X = np.log(normalized_pressure)
    
    f = interpolate.interp1d(x,y,fill_value="extrapolate")
    normalized_temperature = f(X)   
    
    return normalized_temperature

def calc_xsec_data(molecule,component,numin,numax,P,T,step,d_path,r_path,remake=False):
    
    name = "%s/%s_%s_%s_%s_%s_%s.npy"%(r_path,molecule,numin,numax,P,T,step)
    
    if not os.path.isfile(name) or remake:
        print name
        
        calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step,float(P),float(T))
        data = calc.read_data()
      
        
    else:
        nu, coef = np.load(name)
    
    return nu,coef

def load_xsec_data(molecule_list,Tinput,Pinput,name="test",remake=False):
    
    numin = 400
    numax = 10000
    step = 1    
    nu = np.arange(numin,numax,step)
    
    if os.path.isfile("%s.npy"%name) and not remake:
        all_normalized_xsec = np.load("%s.npy"%name)
        return all_normalized_xsec
    
    all_normalized_xsec = []
    
    for molecule in molecule_list:
        try:
            component = [HITRAN_Match[molecule],1,1]
        except:
            all_normalized_xsec.append(np.zeros(shape=(len(Tinput),len(nu))))
            continue
        
        d_path      = os.path.join(HITRAN_Lines,molecule)
        r_path      = "../../input/calculated_xsec/%s"%molecule
        
        if not os.path.isdir(r_path):
            os.makedirs(r_path)
    
        # PArray need to be in increasing order and linear space for interpolation to save time
        PArray = [200000,100000,10000,1000,100,10,1,0.1,0.01,0.001][::-1]
        TArray = [150,175,200,225,250,275,300,325,350]
        
        data = np.zeros(shape=(len(PArray),len(TArray),len(nu)))
        
        for i,P in enumerate(PArray):
            for j,T in enumerate(TArray):
                try:
                    nu, sigma_cm = calc_xsec_data(molecule,component,numin,numax,P,T,step,d_path,r_path)
                except:
                    print "error"
                data[i][j] = sigma_cm  
                
        my_interpolating_function = RegularGridInterpolator((np.log10(PArray),TArray, nu), data)
        
        normalized_xsec = []
        for P_E,T_E in zip(Pinput,Tinput):
            if P_E < 0.001:
                P_E = 0.001
            #    continue
            #print "%.4g"%T_E,"%.4g"%np.log10(P_E)
            if T_E < 150:  # need to calculate for lower temeprature than 175
                T_E = 150

            pts = []
            for n in nu:
                pts.append([np.log10(P_E),T_E,n])
            
            normalized_xsec.append(my_interpolating_function(pts))
        
        all_normalized_xsec.append(normalized_xsec)
    
    np.save("%s.npy"%name, np.array(all_normalized_xsec))
    
    return np.array(all_normalized_xsec)



def load_xsec_data_complex(molecule_list,Tinput,Pinput,wncon,name="test",remake=False):

    wn = []
    for i in wncon:
        numin,numax,step = i
        wn = np.concatenate([wn,np.arange(numin,numax,step)[:-1]])
    
    if os.path.isfile("%s.npy"%name) and not remake:
        all_normalized_xsec = np.load("%s.npy"%name)
        return all_normalized_xsec
    
    all_normalized_xsec = []
    
    for molecule in molecule_list:
        try:
            component = [HITRAN_Match[molecule],1,1]
        except:
            print molecule, "failed"
            all_normalized_xsec.append(np.zeros(shape=(len(Tinput),len(wn))))
            continue
        
        d_path      = os.path.join(HITRAN_Lines,molecule)
        r_path      = "../../input/calculated_xsec2/%s"%molecule
        
        if not os.path.isdir(r_path):
            os.makedirs(r_path)
    
        # PArray need to be in increasing order and linear space for interpolation to save time
        PArray = [200000,100000,10000,1000,100,10,1,0.1,0.01,0.001][::-1]
        TArray = [150,175,200,225,250,275,300,325,350]
        
        data = np.zeros(shape=(len(PArray),len(TArray),len(wn)))
        
        for i,P in enumerate(PArray):
            for j,T in enumerate(TArray):
                try:
                    sigma_cm = []
                    for wnin in wncon:
                        numin,numax,step = wnin
                        nu, sigma_ = calc_xsec_data(molecule,component,numin,numax,P,T,step,d_path,r_path)
                        sigma_cm = np.concatenate([sigma_cm,sigma_[:-1]])
                    
                    
                except:
                    print "error"
                data[i][j] = sigma_cm  
                
        my_interpolating_function = RegularGridInterpolator((np.log10(PArray),TArray, wn), data)
        
        normalized_xsec = []
        for P_E,T_E in zip(Pinput,Tinput):
            if P_E < 0.001:
                P_E = 0.001
            #    continue
            #print "%.4g"%T_E,"%.4g"%np.log10(P_E)
            if T_E < 150:  # need to calculate for lower temeprature than 175
                T_E = 150

            pts = []
            for n in wn:
                pts.append([np.log10(P_E),T_E,n])
            
            normalized_xsec.append(my_interpolating_function(pts))
        
        all_normalized_xsec.append(normalized_xsec)
    
    np.save("%s.npy"%name, np.array(all_normalized_xsec))
    
    return np.array(all_normalized_xsec)




def load_hitran_processed_isoprene():
    
    datapath = "../../input/absorption_data/HITRAN_Cross_Section/isoprene/C5-H8_323.1K-760.0K_600.0-6500.0_0.11_N2_505_43.xsc"
    
    file = open(datapath,"r").read().split("\n")
    
    header = file[0]
    headerinfo = header.split()
    
    numin   = headerinfo[1]
    numax   = headerinfo[2]
    npts    = headerinfo[3]
    T       = headerinfo[4]
    P       = headerinfo[5]
    maxres  = headerinfo[6]
    molecule= headerinfo[7]
    broaden = headerinfo[8]
    note    = headerinfo[9]
    max = maxres[:-5]
    res = maxres[-5:]
    
    ylist = np.array(np.concatenate([i.split() for i in file[1:]]),dtype=float)
    xlist = np.linspace(float(numin),float(numax),len(ylist))
    
    return xlist,ylist

def blackbody_lam(lam, T):
    """ 
    Blackbody as a function of wavelength (um) and temperature (K).
    """
    lam = 1e-6 * lam # convert to metres
    return 2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))

def blackbody_nu(wn,T): # nu is frequency (Hz) , T is temperature (K)
    """ 
    Blackbody as a function of wavenumber (cm^-1) and temperature (K).
    """
    nu = wn*3e10 # convert wavenumber to frequency
    return (2*h*nu**3/c**2)*(1/(np.exp((h*nu)/(k*T))-1))

def blackbody_nu2(nu,T): # nu is frequency (Hz) , T is temperature (K)

    h=6.6244e-34 # Planck constant in Js
    c=2.998e8     # speed of light in m/s
    k=1.38e-23    #Boltzmann constant in J/K%
    
    u_c=(8*np.pi*nu**2*k*T)/c**3 #classical part of energy density
    x=(h*nu)/(k*T)
    eb=x/(np.exp(x)-1 ) #Einstein-Bose correction factor - no units
    u_q=u_c*eb  # energy density
    S=u_q*c/4.0 # radiation
    return S

def test_bb():

    numin = 10
    numax = 10000
    step = 1    
    wn = np.arange(numin,numax,step)     
    
    plt.plot(wn,blackbody_nu(wn,300))
    plt.show()

    
def emission_simulation():
    
    numin = 400
    numax = 10000
    step = 1    
    wn = np.arange(numin,numax,step)    
    
    nurange = wn*(3e10)
    
    
    Z_E,P_E,T_E = read_modtran_input()
    Z_E = Z_E[1:]-Z_E[:-1]
    P_E = (P_E[1:]+P_E[:-1])*0.5
    T_E = (T_E[1:]+T_E[:-1])*0.5
    
    molecule_list,molecule_MR,MR_pressure = load_Earth_MR()
    normalized_molecule_MR = interpolate_MR(molecule_MR,MR_pressure,P_E)    
    
    path = "../../input/temporary/em_xsec"
    #xsec_name = os.path.join(path,"%s_%s_%s"%(main_molecule,stellar,"Earth"))
    xsec_name = os.path.join(path,"Earth_Test4")
    molecule_xsec = load_xsec_data(molecule_list,T_E,P_E,name=xsec_name,remake=True)   
    
    
    count = np.arange(len(P_E))
    
    
    Total_Intensity = np.zeros(len(wn))

    tau_of_each_layer = []
    
    Surface_Intentsity = blackbody_nu(wn,300)    
    
    for C,T,P,L in zip(count,T_E,P_E,Z_E):
        tau_per_layer = np.zeros(len(wn))
        for moi, molecule in enumerate(molecule_list):
            
            sigma_cm = molecule_xsec[moi][C]
            sigma = sigma_cm*0.0001
            n = P/(k*T)*normalized_molecule_MR[moi][C]
            tau_per_molecule = n*sigma*L
            
            tau_per_layer+=tau_per_molecule
            
        
        Surface_Intentsity *= np.exp(-tau_per_layer)    
    
        tau_of_each_layer.append(tau_per_layer)
        
    Total_Intensity+=Surface_Intentsity

    # I need to setup a model that goes from the top of the atmosphere and work my way down
    # at each layer, the 
    for C,T,P,L in zip(count,T_E,P_E,Z_E):

        tau_cur_layer = tau_of_each_layer[C]
        Layer_Intentsity = blackbody_nu(wn,T)*tau_cur_layer
        
        # C+1 because itshould not be attenuated by itself.
        #if C < 20:#count[-1]:
        for j in count[C:]:
            Layer_Intentsity *= np.exp(-tau_of_each_layer[j])
        Total_Intensity += Layer_Intentsity


    #plt.yscale("log")
    plt.plot(wn,blackbody_nu(wn,300))
    plt.plot(wn,blackbody_nu(wn,280))
    plt.plot(wn,blackbody_nu(wn,260))
    plt.plot(wn,blackbody_nu(wn,240))
    plt.plot(wn,blackbody_nu(wn,200))
    plt.plot(wn,blackbody_nu(wn,300)-Total_Intensity)
    #plt.ylim(1e2,1e8)
    #plt.xlim(400,2000)
    plt.show()
    


if __name__ == "__main__":
    emission_simulation()
    #test_bb()
    

