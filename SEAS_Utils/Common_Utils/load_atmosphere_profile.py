"""

This is a submodule for the data loader to load atmosphere profiles
Should not be called by user

need to add a way so that users can create their own input

not sure if helper functions should belong in astrophysics or here

"""
import os
import sys
import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Utils.System_Utils.optimization as opt
import SEAS_Main.Physics.astrophysics as calc


VERBOSE = False

def load_atmosphere_pressure_layers(user_input):
        
    Surface_Pressure = float(user_input["Atmosphere"]["Surface_Pressure"])
    P_Cutoff         = float(user_input["Atmosphere"]["TS_P_Cut_Off"])

    if Surface_Pressure == -1:
        print("undefined surface pressure, set to 100000 Pa")
        Surface_Pressure = 100000

    normalized_pressure = []
    P = Surface_Pressure
    
    while P > P_Cutoff:
        normalized_pressure.append(float("%.3g"%P))
        P = P*np.e**-(1./float(user_input["Atmosphere"]["Sub_Layers"]))
        

    return normalized_pressure

def interpolate_atmosphere_profile(user_input):
    """
    This is a simplistic way of calculating the scale height... 
    need to double check if much differ from the real way of doing this
    may have some trouble for much higher up
    """
    
    TP_Pressure         = user_input["Prototype"]["Input_Pressure"]
    TP_Temperature      = user_input["Prototype"]["Input_Temperature"]
    normalized_pressure = user_input["Prototype"]["Normalized_Pressure"]
    Molecule_List       = user_input["Prototype"]["Molecule_List"] 
    MR_Profile          = user_input["Prototype"]["Input_MR_Profile"] 
    
    assert np.log10(TP_Pressure)[0] < 100 and np.log10(TP_Pressure)[-1] > -100
    
    x = np.concatenate([[100],np.log10(TP_Pressure),[-100]])
    y = np.concatenate([[TP_Temperature[0]],TP_Temperature,[TP_Temperature[-1]]])
    
    
    normalized_MR_profile = MR_Profile.copy()
    for molecule in Molecule_List:
        profile = MR_Profile[molecule]
        k = np.concatenate([[profile[0]],profile,[profile[-1]]])
        normalized_MR_profile[molecule] = interp1d(x,k)(np.log10(normalized_pressure))
    
    user_input["Prototype"]["Normalized_Temperature"]  = interp1d(x,y)(np.log10(normalized_pressure))
    user_input["Prototype"]["Normalized_MR_Profile"]   = normalized_MR_profile

    return user_input

# Individual functions 

def calculate_scale_height(user_input):
    
    Molecule_List          = user_input["Prototype"]["Molecule_List"] 
    normalized_MR_profile  = user_input["Prototype"]["Normalized_MR_Profile"]
    normalized_temperature = user_input["Prototype"]["Normalized_Temperature"]
    Surface_Gravity        = user_input["Planet"]["Surface_Gravity"]
    
    molecular_weight_list = calc.get_MolWeight(Molecule_List)
    profile_per_layer = np.array(normalized_MR_profile.values()).T
   
    normalized_scale_height = np.zeros(len(normalized_temperature))
    normalized_mean_mw = np.zeros(len(normalized_temperature))
    for i,abundance_per_layer in enumerate(profile_per_layer):
        mean_mw      = calc.calc_MeanMolWeight(abundance_per_layer, molecular_weight_list)
        scale_height = calc.calc_H(normalized_temperature[i], mean_mw*mH, Surface_Gravity)
        normalized_scale_height[i] = scale_height
        normalized_mean_mw[i] = mean_mw
      
    user_input["Prototype"]["Normalized_Scale_Height"] = normalized_scale_height
    user_input["Prototype"]["normalized_mean_mw"] = normalized_mean_mw
    
    return user_input


class Atmosphere_Profile_Loader():
    
    def __init__(self):
        """ will change when user_input loading is updated"""
        pass

    def load_Atmosphere_Profile_from_Photochemistry_Code(self,user_input,Profile):
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
        
        output_data=np.genfromtxt(Profile, skip_header=2, skip_footer=0, unpack=False, encoding='ascii') 
    
        z_centers=output_data[:,0]*km2m # Center of altitude bins, km converted to cm
        z_lower=output_data[:,1]*km2m # Lower edge of altitude bins, km converted to cm
        z_upper=output_data[:,2]*km2m # Upper edge of altitude bins, km converted to cm
        T_z = output_data[:,3] # Temperature(z), in K
        P_z = output_data[:,4] # Pressure(z), in Pa 
        n_z_species=output_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)
        n_z=np.sum(n_z_species,1) #sum number densities across species. This is a profile for the whole atmosphere.
            
        # Truncating the TP profile by imposing a pressure cutoff
        # The flat spectra problem could be due to the fact that the lowest pressure is 0.03 Pa
        category  = user_input["Spectra"]["Category"]
        threshold = float(user_input["Atmosphere"]["%s_P_Cut_Off"%category])
        accepted  = len([x for x in P_z if x > threshold])
        P_z = P_z[:accepted]
        T_z = T_z[:accepted]
        z_centers = z_centers[:accepted]
            
        MR_Profile = {}
        Molecule_List = []
        for i,info in enumerate(n_z_species.T):
            l = key_num2species[i+1]
            if max(info/n_z) > float(user_input["Prototype"]["Threshold"]): #default is 10**-7
                MR_Profile[str(l)] = (info/n_z)[:len(P_z)]
                Molecule_List.append(l)
    
        user_input["Prototype"]["Molecule_List"]      = Molecule_List
        user_input["Prototype"]["Input_MR_Profile"]   = MR_Profile
        user_input["Prototype"]["Input_Pressure"]     = P_z
        user_input["Prototype"]["Input_Temperature"]  = T_z
        
        
        return user_input
    
    def load_Atmosphere_Profile_from_Boxcar(self,user_input,Profile):
        
        user_input["Prototype"]["Molecule_List"]           = Profile["Molecule"].keys()
        user_input["Prototype"]["Normalized_MR_Profile"]   = Profile["Molecule"]
        user_input["Prototype"]["Normalized_Pressure"]     = Profile["Pressure"]
        user_input["Prototype"]["Normalized_Temperature"]  = Profile["Temperature"]
        user_input["Prototype"]["Normalized_Scale_Height"] = Profile["Scale_Height"]
        
        return user_input
    
    def load_Atmosphere_Profile_from_CCM(self,user_input,Profile):
        
        Molecule_List = user_input["Prototype"]["Source_Header"]
        output_data   = np.genfromtxt(Profile, unpack=False, encoding='ascii') 
    
        P_z = output_data[:,0][::-1]*100  # hPa to Pa, also reverse the order so that it is from bottom up
        T_z = output_data[:,1][::-1]
        n_z_species = np.flip(output_data[:,2:],0)
        n_z=np.sum(n_z_species,1)
        
        # Truncating the TP profile by imposing a pressure cutoff
        # The flat spectra problem could be due to the fact that the lowest pressure is 0.03 Pa
        category  = user_input["Spectra"]["Category"]
        threshold = float(user_input["Atmosphere"]["%s_P_Cut_Off"%category])
        accepted  = len([x for x in P_z if x > threshold])
        P_z = P_z[:accepted]
        T_z = T_z[:accepted]
        
        # This is a hacky thing that needs to be changed
        T_z = [x if x < 400 else 400 for x in T_z]
    
        MR_Profile = {}
        for molecule,data in zip(Molecule_List,n_z_species.T):
            MR_Profile[molecule] = (data/n_z)[:len(P_z)]
    
        user_input["Prototype"]["Molecule_List"]      = Molecule_List
        user_input["Prototype"]["Input_MR_Profile"]   = MR_Profile
        user_input["Prototype"]["Input_Pressure"]     = P_z
        user_input["Prototype"]["Input_Temperature"]  = T_z
        
        return user_input
    
    def load_Atmosphere_Profile_from_Earth(self,user_input):
        """
        no need to specify scenario file in this case
        """
        
        MR_File = "../../SEAS_Input/Atmosphere_Data/MR_Profile/earth.txt"
        TP_File = "../../SEAS_Input/Atmosphere_Data/TP_profile/earth.txt"
        
        data = multi_column_file_loader(MR_File,type="mixed")
        
        Molecule_List = []
        MR_Profile    = {}
        for i in data[1:]:
            molecule = i[0]
            Molecule_List.append(molecule)
            MR_Profile[molecule] = [float(x)/100. for x in i[1:]]
    
        MR_pressure = [float(x) for x in data[0][1:]]
        
        P,T = np.genfromtxt(TP_File).T
     
        MR_temperature = interp1d(np.log10(P),T)(np.log10(MR_pressure))
     
     
        user_input["Prototype"]["Molecule_List"]      = Molecule_List
        user_input["Prototype"]["Input_MR_Profile"]   = MR_Profile
        user_input["Prototype"]["Input_Pressure"]     = MR_pressure
        user_input["Prototype"]["Input_Temperature"]  = MR_temperature
        
    
        return user_input
    
