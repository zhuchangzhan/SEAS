import os
import sys
import h5py
import hashlib
import numpy as np
from numba import jit
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Main.Physics.astrophysics as calc
from SEAS_Main.Simulation.cloud import Simple_Gray_Cloud_Simulator
import SEAS_Utils.System_Utils.optimization as opt
import SEAS_Utils.Common_Utils.db_management2 as dbm
import SEAS_Utils.Common_Utils.interpolation as interp

@opt.timeit
def load_Observation_Data(user_input):
    return

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
    
    user_input["Planet"]["Surface_Gravity"] = calc.calc_SurfaceG(M_planet, R_planet)
    user_input["Atmosphere"]["Base_TS_Value"]   = (R_planet/R_Star)**2
    
    return user_input


#@jit(nopython=True)

@opt.timeit
def load_Atmosphere_Profile(user_input,source,scenario_file=None):
    """
    Toggler for loading atmosphere profile, which includes:
    1. TP Profile
    2. MR Profile
    
    By Default is loading source from photochemistry code:
    """
    if source == "Photochemistry" and scenario_file != None:
        
        # Create a Hash of the TP profile so that cross sections only need to be generated once per profile
        # Some testing here or a better format is needed
        return load_Atmosphere_Profile_from_Photochemistry_Code(user_input,scenario_file)
    else:
        print("Not Implemented")
        sys.exit()

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
    for j,info in enumerate(n_z_species.T):
        l = key_num2species[j+1]
        if max(info/n_z) > float(user_input["Prototype"]["Threshold"]): #default is 10**-7
            MR_Profile[str(l)] = (info/n_z)[:len(P_z)]
            Molecule_List.append(l)

    user_input["Prototype"]["Molecule_List"] = Molecule_List
    user_input["Prototype"]["MR_Profile"] = MR_Profile
    user_input["Prototype"]["Normalized_Scale_Height"] = z_centers
    user_input["Prototype"]["TP_Profile"]["Normalized_Pressure"] = P_z
    user_input["Prototype"]["TP_Profile"]["Normalized_Temperature"] = T_z
    
    return user_input

@opt.timeit
def load_Absorption_Cross_Section(user_input,reuse=True):
    
    # Load Absorption cross section
    info = Xsec_Loader(user_input,reuse)
    for molecule in user_input["Prototype"]["Molecule_List"]:
        info.load_HITRAN(molecule)
        user_input["Xsec"]["Molecule"][molecule] = info.xsec[molecule]
    
    # Load Rayleigh Scattering cross section
    user_input["Xsec"]["Rayleigh"]["Value"] = info.load_rayleigh_scattering(user_input["Prototype"]["Molecule_List"])
    user_input["Xsec"]["Cloud"]["Value"]    = info.load_gray_cloud()
    user_input["Xsec"]["nu"]                = info.nu
    return user_input


class Xsec_Loader():
    """
    Absorption Cross Section Loader.
    
    is there a way to get rid of needing in have user_input
    loaded in "load_HITRAN"? it's already loaded with __init__
    """

    def __init__(self, user_input,reuse=True):
        self.user_input = user_input
        self.reuse = reuse
        
        TP_Profile = user_input["Prototype"]["TP_Profile"]
        self.normalized_pressure        = TP_Profile["Normalized_Pressure"]
        self.normalized_temperature     = TP_Profile["Normalized_Temperature"]
        self.normalized_scale_height    = user_input["Prototype"]["Normalized_Scale_Height"]
        self.xsec = {}
        
        self.DB_DIR = "../../SEAS_Input/Cross_Section/HDF5_DB"
    
        # for interpolating the cross section
        self.x = self.user_input["Xsec"]["Molecule"]["T_Grid"]
        self.X = self.normalized_temperature 
        
    def load_wavelength(self):
        pass
    
    def load_HITRAN(self, molecule, savepath=None, savename=None):
        # This step can be optimized for retrieval by only loading it once
    
        hash = self.user_input["Data_IO"]["Hash"]
        savepath = savepath or "../../SEAS_Input/Cross_Section/Generated/%s"%hash
        savename = savename or "%s_%s.npy"%(molecule,hash)
        filepath = os.path.join(savepath,savename)
    
        if self.reuse and os.path.isfile(filepath):
            self.nu,self.xsec[molecule] = np.load(filepath)
            print("%s Cross Section Loaded"%molecule)
        else:
            nu = h5py.File("%s/%s.hdf5"%(self.DB_DIR,"nu"), "r")
            self.nu = np.array(nu["results"])
            
            xsec = h5py.File("%s/%s.hdf5"%(self.DB_DIR,molecule), "r")
            raw_cross_section_grid = np.array(xsec["results"])
             
            self.xsec[molecule] = self.grid_interpolate(raw_cross_section_grid)
            np.save(filepath, [self.nu,self.xsec[molecule]])
            print("%s Cross Section Saved"%molecule)

    def load_Exomol(self, molecule):
        
        self.xsec = ["Exomol_%s not implemented"%molecule]
            
    def load_NIST(self, molecule):
        
        self.xsec = ["NIST_%s not implemented"%molecule] 

    @opt.timeit
    def grid_interpolate(self,xsec_grid):
        
        T_Grid = np.array(self.user_input["Xsec"]["Molecule"]["T_Grid"],dtype=float)
        P_Grid = np.array(self.user_input["Xsec"]["Molecule"]["P_Grid"],dtype=float)
        
        # creating the 2D interpolation function
        # the [::-1] is because it needs to be strictly ascending
        f = RegularGridInterpolator((np.log10(P_Grid[::-1]),T_Grid, self.nu),xsec_grid[::-1])
        wave = len(self.nu)
        layer = len(self.normalized_pressure)
        normalized_xsec = np.zeros((layer,wave))
        for i,(P_E,T_E) in enumerate(zip(np.log10(self.normalized_pressure),self.normalized_temperature)):
            normalized_xsec[i] = f(np.array([np.ones(wave)*P_E, np.ones(wave)*T_E, self.nu]).T)
            
        return normalized_xsec

    def load_rayleigh_scattering(self,molecules):
        """
        currently not caring about biosignature molecule rayleigh?
        """
        
        Rayleigh_array = {}
        for molecule in molecules:
            Rayleigh_array[molecule] = calc.calc_rayleigh(molecule, self.nu)
        print("Rayleigh Scatter Loaded")
        return Rayleigh_array

    def load_gray_cloud(self):

        normalized_cloud_xsec = []
        cloud_deck    = float(self.user_input["Xsec"]["Cloud"]["Deck"])
        cloud_opacity = float(self.user_input["Xsec"]["Cloud"]["Opacity"])
        normalized_pressure = np.array(self.user_input["Prototype"]["TP_Profile"]["Normalized_Pressure"],dtype=float)
        
        #c = Simple_Gray_Cloud_Simulator(cloud_deck,cloud_opacity)
        
        for i,P in enumerate(normalized_pressure):
            if P < cloud_deck:
                normalized_cloud_xsec.append(np.zeros(len(self.nu)))
            else:
                normalized_cloud_xsec.append(np.ones(len(self.nu))*cloud_opacity)
        print("Gray Cloud Loaded")
        return normalized_cloud_xsec   




















    
    

