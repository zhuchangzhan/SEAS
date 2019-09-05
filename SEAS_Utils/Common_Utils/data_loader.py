
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *


def load_Atmosphere_Profile(user_input,source,scenario_file=None):
    """
    Toggler for loading atmosphere profile, which includes:
    1. TP Profile
    2. MR Profile
    
    By Default is loading source from photochemistry code:
    """
    if source == "Photochemistry" and scenario_file != None:
        return load_Atmosphere_Profile_from_Photochemistry_Code(user_input,scenario_file)

def load_Atmosphere_Profile_from_Photochemistry_Code(user_input,scenario_file):
    """
    generate mixing ratio profile from renyu's simulation
    
    stellar_type = which type of star to consider: Sun, Ma, Mq
    main_gas = dominant composition of the atmosphere
    base_mr = base mixing ratio of biosignature gas in consideration
    mu_atm = average molar mass of atmosphere for each scenario (grams)
    g_atm = gravitational constant of atmosphere for each scenario (cm s**-2)
    """
    ########################
    ###Read in key file, to translate between Hu code numbering convention and English names.
    ########################
    species_name_file= user_input["Data_IO"]["File_Path"]["Species_Name"] #Directory storing the key mapping species to the integers 1-111 used to ID them in the code
    imported_key = np.genfromtxt(species_name_file, skip_header=1, skip_footer=0, delimiter='\t', dtype=None, encoding='ascii') #Import mapping between numerical ID in code and species name.
     
    key_species2num={} #initalize empty dictionary
    key_num2species={} #initalize empty dictionary, to hold mapping FROM species name TO species number
    
    for i in range(0, len(imported_key)):
        species_number, species_name = imported_key[i]
        key_species2num[species_name]=species_number #Holds mapping FROM species name TO species number
        key_num2species[species_number]=species_name #Holds mapping FROM species number TO species name
    
    output_data=np.genfromtxt(scenario_file, skip_header=2, skip_footer=0, unpack=False, encoding='ascii') 

    z_centers=output_data[:,0]*km2cm # Center of altitude bins, km converted to cm
    z_lower=output_data[:,1]*km2cm # Lower edge of altitude bins, km converted to cm
    z_upper=output_data[:,2]*km2cm # Upper edge of altitude bins, km converted to cm
    T_z=output_data[:,3] # Temperature(z), in K
    P_z=output_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    P_z_Pa = output_data[:,4] # Pressure(z), in Pa 
    n_z_species=output_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)
    n_z=np.sum(n_z_species,1) #sum number densities across species. This is a profile for the whole atmosphere.
        
    # Saving the TP profile
    # This code can be optimized but not urgent
    newP,newT,z_out = [],[],[]
    category = user_input["Spectra"]["Category"]
    for P,T,Z in zip(P_z_Pa,T_z,z_centers/100): # zcenter in meter
        if P < float(user_input["Atmosphere"]["%s_P_Cut_Off"%category]):
            break
        newT.append(T)
        newP.append(P)
        z_out.append(Z)
    TP_Profile = [newP,newT,z_out]
    
    MR_Profile = {}
    Molecule_List = []
    for j,info in enumerate(n_z_species.T):
        l = key_num2species[j+1]
        Molecule_List.append(l)
        if max(info/n_z) > float(user_input["Prototype"]["Threshold"]): #default is 10**-7
            MR_Profile[str(l)] = (info/n_z)[:len(newP)]

    user_input["Prototype"]["Molecule_List"] = Molecule_List
    user_input["Prototype"]["TP_Profile"] = TP_Profile
    user_input["Prototype"]["MR_Profile"] = MR_Profile
    
    return user_input





