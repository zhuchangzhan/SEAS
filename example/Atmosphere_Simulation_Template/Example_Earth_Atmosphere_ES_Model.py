"""

Simulated Earth Emission Spectra Model

validated using MODTRAN 



"""


import os
import sys
import tqdm
import glob
import hashlib
from scipy import stats
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Utils.Common_Utils.configurable as config
import SEAS_Utils.System_Utils.optimization as opt
import SEAS_Utils.Common_Utils.data_loader as load 
import SEAS_Main.Physics.astrophysics as calc 
import SEAS_Main.Simulation.emission_spectra_simulator as ES
from SEAS_Main.Physics.noise import Photon_Noise


user_input = config.Configuration("../../config/user_input.cfg")
VERBOSE = bool(user_input["Data_IO"]["Logging"]["VERBOSE"])

@opt.timeit
def Generate_Atmosphere_Spectra(user_input):
    
    simulation = ES.Emission_Spectra_Simulator(user_input)
    simulation.load_atmosphere_geometry_model()
        
    return simulation.user_input

@opt.timeit
def Simulate_Atmosphere_Observation(user_input):
       
    nu_               = user_input["Spectra"]["Wavelength"]
    Total_Intensity   = user_input["Spectra"]["Total_Intensity"]
    
    noise = Photon_Noise(user_input)
    bin_edges, bin_width, bin_centers = noise.determine_bin()
    nu,  convolved_Height = noise.calculate_convolve(nu_, Total_Intensity)    
    Bin_Intensity, bin_edges, binnumber = stats.binned_statistic(10000./nu_[::-1], Total_Intensity[::-1], bins=bin_edges)

    T_Star = 4000
    R_Star = 0.6*R_Sun
    R_Planet = 1*R_Earth
    Planet_Star_A_Ratio = R_Planet**2/R_Star**2
    Star_BB_Bin         = calc.blackbody_nu(10000./bin_centers,T_Star)

    user_input["Spectra"]["bin_width"]          = bin_width
    user_input["Spectra"]["bin_centers"]        = bin_centers
    user_input["Spectra"]["nu"]                 = nu
    user_input["Spectra"]["bin_values"]         = convolved_Height
    user_input["Spectra"]["bin_intensity"]      = Bin_Intensity
    user_input["Spectra"]["Planet_Star_Ratio"]  = Bin_Intensity/Star_BB_Bin*Planet_Star_A_Ratio

    return user_input

    """
    main_mol    = user_input["Prototype"]["Atmosphere_Type"]
    stellar     = user_input["Prototype"]["Stellar_Type"]
    R_Planet    = float(user_input["Planet"]["R_Planet"])*R_Earth
    R_Star      = float(user_input["Star"]["R_Star"])*R_Sun
    T_Star      = float(user_input["Star"]["T_Star"])
    Distance    = float(user_input["System"]["D_Star_Observer"])*Psec  
    
    D_obs       = float(user_input["Telescope"]["Aperture"])
    R_obs       = D_obs/2
    Duration    = float(user_input["Telescope"]["Duration"])*3600. # hr to second
    Quantum     = float(user_input["Telescope"]["Quantum_Efficiency"])
    Noise_M     = float(user_input["Telescope"]["Noise"]["multiplier"])
    
    # calculate number of photons 
    B_Body      = calc.blackbody_lam(bin_centers*10**-6, T_Star)
    Bin_width   = bin_width*10**-6
    A_Star      = np.pi*R_Star**2
    Psi_Tele    = np.pi*R_obs**2/Distance**2
    E_Total     = B_Body*Bin_width*A_Star*Psi_Tele*Duration
    star_photon = (E_Total*bin_centers*10**-6)/(HPlanck*CLight)*Quantum
    
    user_input["Spectra"]["Signal"] = R_Planet**2/R_Star**2*star_photon
    user_input["Spectra"]["Noise"]  = Noise_M*np.sqrt(star_photon)
    user_input["Spectra"]["SNR"]    = user_input["Spectra"]["Signal"]/user_input["Spectra"]["Noise"]
    """
    return user_input
    

def Emission_Forward_Model_Architecture():
    
    user_input = config.Configuration("../../config/user_input.cfg")
    user_input["Prototype"]["Source"]        = "Earth"
    user_input["Prototype"]["Source_Header"] = "H2O CO2 CH4 O2 O3 N2".split(" ")
    user_input["Config"]["molecule_turnoff"] = None
    
    file = glob.glob("../../SEAS_Input/Atmosphere_Data/MR_Profile/*.txt")[0]
    
    
    # Load Surface_Gravity, Base_TS_Value
    user_input = load.load_Astrophysical_Properties(user_input)
    
    # Loading TP_Profile, MR_Profile, Molecule List, 
    user_input = load.load_Atmosphere_Profile(user_input, scenario_file=file)
 
    # Load absorption cross section for all molecule and effect (CIA, Cloud, etc)
    user_input = load.load_Absorption_Cross_Section(user_input,True)
    
    # Load Atmosphere model and generate theoretical spectra
    user_input = Generate_Atmosphere_Spectra(user_input)
    
    # Load Observation parameters and generate simulated observation
    user_input = Simulate_Atmosphere_Observation(user_input)

    #nu = user_input["Spectra"]["nu"] 
    a1 = user_input["Spectra"]["bin_centers"]    
    b1 = user_input["Spectra"]["bin_intensity"]      
    bb = user_input["Spectra"]["Planet_Star_Ratio"]  
    
    nu = 10000./a1
    
    plt.plot(nu, calc.blackbody_nu(nu,320),color="k")
    plt.plot(nu, calc.blackbody_nu(nu,300),color="0.9")
    plt.plot(nu, calc.blackbody_nu(nu,280),color="0.7")
    plt.plot(nu, calc.blackbody_nu(nu,260),color="0.6")
    plt.plot(nu, calc.blackbody_nu(nu,240),color="0.4")
    plt.plot(nu, calc.blackbody_nu(nu,220),color="0.2")
          
    plt.plot(nu, b1)
        
    plt.gca().invert_xaxis()
    plt.xlim(400,1600)
    plt.show()

    

if __name__ == "__main__":
    
    Emission_Forward_Model_Architecture()