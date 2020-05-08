"""
Here is a construction site for building SEAS.
Will have to clean up once finished
Let's build the user input section first
"""

import os
import sys
import tqdm
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
import SEAS_Main.Simulation.transmission_spectra_simulator as TS
from SEAS_Main.Physics.noise import Photon_Noise


user_input = config.Configuration("../../config/user_input.cfg")
VERBOSE = bool(user_input["Data_IO"]["Logging"]["VERBOSE"])

@opt.timeit
def Generate_Atmosphere_Spectra(user_input):
    
    simulation = TS.Transmission_Spectra_Simulator(user_input)
    simulation.load_boxcar_model()
    return simulation.user_input


@opt.timeit
def Forward_Boxcar_Model_Architecture():
    
    global user_input
    
    user_input["Prototype"]["Source"]="Boxcar"
    
    # use percent for consistency for now
    Profile = {}
    Profile["Molecule"] = {}
    
    # even though it's only 1 layer, still use list to keep code consistent
    Profile["Temperature"]      = [300]
    Profile["Pressure"]         = [100000]
    Profile["Molecule"]["H2O"]  = [1]
    Profile["Molecule"]["CH4"]  = [1]
    Profile["Molecule"]["N2"]   = [100 - Profile["Molecule"]["H2O"][0] - Profile["Molecule"]["CH4"][0]]
    
    
    # Loading TP_Profile, MR_Profile, Molecule List, 
    user_input = load.load_Atmosphere_Profile(user_input, scenario_file=Profile)
    
    # Load absorption cross section for all molecule and effect (CIA, Cloud, etc)
    user_input = load.load_Absorption_Cross_Section(user_input,False)
    
    # Load Atmosphere model and generate theoretical spectra
    user_input = Generate_Atmosphere_Spectra(user_input)
    
    plt.plot(10000./user_input["Spectra"]["Wavelength"],
             user_input["Spectra"]["Total_Transit_Signal"])
    plt.xscale("log")
    plt.show()
    

if __name__ == "__main__":
    
    
    
    Forward_Boxcar_Model_Architecture()
















