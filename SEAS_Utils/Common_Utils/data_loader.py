"""

Submodule to load all input data for SEAS
serve as main controller 

We need to add checking the molecule availability database



"""

import os
import sys
import hashlib
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Main.Physics.astrophysics as calc
import SEAS_Utils.System_Utils.optimization as opt
import SEAS_Utils.Common_Utils.load_atmosphere_profile as atm_pro
import SEAS_Utils.Common_Utils.load_cross_section as xsec

VERBOSE = False


def multi_column_file_loader(path,spliter=None,type="float",skip=0):
    """
    load data from files that contains multiple columns
    """
    
    with open(path) as data_file:   
        data = [x.split(spliter) for x in data_file.read().split("\n")[skip:]]
        if data[-1] == [] or data[-1] == None or data[-1] == "":
            data = data[:-1]
        
        
        data = [list(x) for x  in zip(*data)]    
        
        if type == "float":
            return np.array(data,dtype=np.float)
        elif type == "int":
            return np.array(data,dtype=np.int)
        elif type == "mixed":
            return data

@opt.timeit
def load_Astrophysical_Properties(user_input):
    """
    There is going to be more development here in the future
    For now it's just loading the star and planet radius, then calculate the baseline
    for the simulation
    """
    
    R_Star          = float(user_input["Star"]["R_Star"])*R_Sun
    R_planet        = float(user_input["Planet"]["R_Planet"])*R_Earth
    M_planet        = float(user_input["Planet"]["M_Planet"])*M_Earth
    

    if user_input["Planet"]["Surface_Gravity"] == "-1":
        user_input["Planet"]["Surface_Gravity"]     = calc.calc_SurfaceG(M_planet, R_planet)
        
    user_input["Atmosphere"]["Base_TS_Value"]   = (R_planet/R_Star)**2
    
    return user_input

@opt.timeit
def load_Atmosphere_Profile(user_input,scenario_file=None):
    """
    Toggler for loading atmosphere profile, which includes:
    1. TP Profile
    2. MR Profile
    
    Example code is loading source from photochemistry code:
    """
    if scenario_file == None:
        print("No Scenario File Specified")
        sys.exit()
        
    info = atm_pro.Atmosphere_Profile_Loader()
    
    source = user_input["Prototype"]["Source"]
    
    
    if source == "Photochemistry":
        # Create a Hash of the TP profile so that cross sections only need to be generated once per profile
        # Some testing here or a better format is needed
        user_input = info.load_Atmosphere_Profile_from_Photochemistry_Code(user_input,scenario_file)
    elif source == "Boxcar":
        user_input = info.load_Atmosphere_Profile_from_Boxcar(user_input,scenario_file)
    elif source == "Earth":
        user_input = info.load_Atmosphere_Profile_from_Earth(user_input)
    elif source == "CCM":
        user_input = info.load_Atmosphere_Profile_from_CCM(user_input,scenario_file)
    else:
        print("Atm. Pro. Method Not Implemented")
        sys.exit()
    
    if source != "Boxcar":
        # SEAS defined its own simulation layers that is different from photochemistry output
        user_input["Prototype"]["Normalized_Pressure"] = atm_pro.load_atmosphere_pressure_layers(user_input)
        user_input = atm_pro.interpolate_atmosphere_profile(user_input)
        user_input = atm_pro.calculate_scale_height(user_input)
        
    return user_input

@opt.timeit
def load_Absorption_Cross_Section(user_input,reuse=True):
    """
    This will get renamed to load_Cross_Section since it's more general than just loading molecular cross section
    
    
    need to rework the molecule cross section loading sequence. 
        create a pandas database which contains data availability
    
    """
    info = xsec.Cross_Section_Loader(user_input,reuse)
    user_input["Xsec"]["nu"]                = info.nu  
    
    
    if user_input["Prototype"]["Source"] in ["Photochemistry","CCM","Earth"]:
        # Hash created for identifying the xsec used for the TP profile
        # this is different for every user
        a = str(user_input["Prototype"]["Normalized_Pressure"])
        b = str(user_input["Prototype"]["Normalized_Temperature"])
        
        if "New" in user_input["Data_IO"]["File_Path"]["DB_DIR"]:
            user_input["Data_IO"]["Hash"] = hashlib.sha224((a+b).encode()).hexdigest()[:9]
        else:
            user_input["Data_IO"]["Hash"] = hashlib.sha224((a+b).encode()).hexdigest()[:8]
        
        # Load Absorption cross section
        for molecule in user_input["Prototype"]["Molecule_List"]:
            info.load_HITRAN(molecule)
            user_input["Xsec"]["Molecule"][molecule] = info.xsec[molecule]
        
        # Load Rayleigh Scattering cross section
        user_input["Xsec"]["Rayleigh"]["Value"] = info.load_rayleigh_scattering(user_input["Prototype"]["Molecule_List"])
        
        # Load Cloud xsec
        # This should always be loaded after the molecular xsec because it needs nu
        if user_input["Xsec"]["Cloud"]["type"] == "grey":
            user_input["Xsec"]["Cloud"]["Value"] = info.load_gray_cloud()
        elif user_input["Xsec"]["Cloud"]["type"] == "Mie":
            user_input["Xsec"]["Cloud"]["Value"] = info.load_mie_cloud()
      
        # Load Collision Induced Absorption
        if "H2" in user_input["Prototype"]["Molecule_List"]:
            user_input["Xsec"]["CIA"]["Enable"] = "True"
            user_input["Xsec"]["CIA"]["H2-H2"]  = info.load_CIA(molecule="H2-H2")
            
            
    elif user_input["Prototype"]["Source"] == "Boxcar":
        
        user_input["Data_IO"]["Hash"] = hashlib.sha224(("Boxcar_300_100000").encode()).hexdigest()[:8] # This is temporary
        
        for molecule in user_input["Prototype"]["Molecule_List"]:
            info.load_HITRAN(molecule)
            user_input["Xsec"]["Molecule"][molecule] = info.xsec[molecule]
    else:
        print("Xsec. Method Not Implemented")
        sys.exit()
                    
    return user_input


    
    


