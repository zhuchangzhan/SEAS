import os
import sys
import h5py
import hashlib
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import miepython as mp
from scipy.interpolate import interp1d, RegularGridInterpolator

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Main.Physics.astrophysics as calc
from SEAS_Main.Simulation.cloud import Simple_Gray_Cloud_Simulator
import SEAS_Utils.System_Utils.optimization as opt
import SEAS_Utils.Common_Utils.db_management2 as dbm
import SEAS_Utils.Common_Utils.interpolation as interp

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
    

    if user_input["Planet"]["Surface_Gravity"] == "-1":
        user_input["Planet"]["Surface_Gravity"]     = calc.calc_SurfaceG(M_planet, R_planet)
        
    user_input["Atmosphere"]["Base_TS_Value"]   = (R_planet/R_Star)**2
    
    return user_input

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
    
def load_Atmosphere_Profile_from_Photochemistry_Code(user_input,Profile):
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

def load_Atmosphere_Profile_from_Boxcar(user_input,Profile):
    
    user_input["Prototype"]["Molecule_List"]           = Profile["Molecule"].keys()
    user_input["Prototype"]["Normalized_MR_Profile"]   = Profile["Molecule"]
    user_input["Prototype"]["Normalized_Pressure"]     = Profile["Pressure"]
    user_input["Prototype"]["Normalized_Temperature"]  = Profile["Temperature"]
    
    return user_input

def load_Atmosphere_Profile_from_CCM(user_input,Profile):
    
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

def load_Atmosphere_Profile_from_Earth(user_input):
    
    
    
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

@opt.timeit
def load_Atmosphere_Profile(user_input,scenario_file=None):
    """
    Toggler for loading atmosphere profile, which includes:
    1. TP Profile
    2. MR Profile
    
    By Default is loading source from photochemistry code:
    """
    
    source = user_input["Prototype"]["Source"]
    
    if scenario_file == None:
        print("Not Implemented")
        sys.exit()
        
    if source == "Photochemistry":
        # Create a Hash of the TP profile so that cross sections only need to be generated once per profile
        # Some testing here or a better format is needed
        user_input = load_Atmosphere_Profile_from_Photochemistry_Code(user_input,scenario_file)
        user_input["Prototype"]["Normalized_Pressure"] = load_atmosphere_pressure_layers(user_input)
        user_input = interpolate_atmosphere_profile(user_input)
        user_input = calculate_scale_height(user_input)
        
    elif source == "Boxcar":
        user_input = load_Atmosphere_Profile_from_Boxcar(user_input,scenario_file)
        
    elif source == "Earth":
        user_input = load_Atmosphere_Profile_from_Earth(user_input)
        user_input["Prototype"]["Normalized_Pressure"] = load_atmosphere_pressure_layers(user_input)
        user_input = interpolate_atmosphere_profile(user_input)
        user_input = calculate_scale_height(user_input)

    elif source == "CCM":
        user_input = load_Atmosphere_Profile_from_CCM(user_input,scenario_file)
        user_input["Prototype"]["Normalized_Pressure"] = load_atmosphere_pressure_layers(user_input)
        user_input = interpolate_atmosphere_profile(user_input)
        user_input = calculate_scale_height(user_input)
    
    else:
        print("Not Implemented")
        sys.exit()
        
    # SEAS defined its own simulation layers that is different from photochemistry output
    # This is to standarize simulation?
    
    return user_input

@opt.timeit
def load_Absorption_Cross_Section(user_input,reuse=True):


    info = Xsec_Loader(user_input,reuse)
    user_input["Xsec"]["nu"]                = info.nu  
    
    
    if user_input["Prototype"]["Source"] in ["Photochemistry","CCM","Earth"]:
        # Hash created for identifying the xsec used for the TP profile
        # this is different for every user
        a = str(user_input["Prototype"]["Normalized_Pressure"])
        b = str(user_input["Prototype"]["Normalized_Temperature"])
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
            user_input["Xsec"]["Cloud"]["Value"]    = info.load_gray_cloud()
        elif user_input["Xsec"]["Cloud"]["type"] == "Mie":
            user_input["Xsec"]["Cloud"]["Value"]    = info.load_mie_cloud()
      
        
    elif user_input["Prototype"]["Source"] == "Boxcar":
        
        for molecule in user_input["Prototype"]["Molecule_List"]:
            print(molecule)
            info.load_HITRAN_single(molecule)
            user_input["Xsec"]["Molecule"][molecule] = info.xsec[molecule]
         
                
        
        
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
        
        
        self.normalized_pressure        = user_input["Prototype"]["Normalized_Pressure"]
        self.normalized_temperature     = user_input["Prototype"]["Normalized_Temperature"]
        
        
        self.DB_DIR = "../../SEAS_Input/Cross_Section/HDF5_DB"
        
        nu = h5py.File("%s/%s.hdf5"%(self.DB_DIR,"nu"), "r")
        self.nu = np.array(nu["results"])
        self.xsec = {}
        
        # for interpolating the cross section
        # self.x = self.user_input["Xsec"]["Molecule"]["T_Grid"]
        # self.X = self.normalized_temperature 
        
    def load_wavelength(self):
        
        self.wave = 10000./self.nu
        
        return self.wave
    
    def load_HITRAN(self, molecule, savepath=None, savename=None,cache=None):
        """
        Loading molecule cross section from database or precalculated
        interpolate cross section along pressure temperature profile
        # This step can be optimized for retrieval by only loading it once if TP profile doesn't change
        """
        
        if cache != None:
            self.xsec[molecule] = cache
    
        hash = self.user_input["Data_IO"]["Hash"]
        savepath = savepath or "../../SEAS_Input/Cross_Section/Generated/%s"%hash
        savename = savename or "%s_%s.npy"%(molecule,hash)
        filepath = os.path.join(savepath,savename)
        
        if self.reuse and os.path.isfile(filepath):
            self.nu,self.xsec[molecule] = np.load(filepath) #self.nu doesn't need to be loaded again
            if VERBOSE:
                print("%s Cross Section Loaded"%molecule)
        else:
            if True:
                print(molecule)
            #try: # This is a temporary solution to add isoprene
                xsec = h5py.File("%s/%s.hdf5"%(self.DB_DIR,molecule), "r")
                raw_cross_section_grid = np.array(xsec["results"])
                self.xsec[molecule] = self.grid_interpolate(raw_cross_section_grid)
            else:
            #except:
                # need to package and modularize this
                # Also change the name of load hitran to load xsec
                xsec = self.load_PNNL(molecule)
                wave = len(self.nu)
                layer = len(self.normalized_pressure)
                normalized_xsec = np.zeros((layer,wave))
                for i in range(layer):
                    normalized_xsec[i] = xsec
                self.xsec[molecule] = normalized_xsec
                print("loaded %s from NIST"%molecule)
            
            if not os.path.isdir(savepath):
                os.makedirs(savepath)   
            np.save(filepath, [self.nu,self.xsec[molecule]]) # no need to save self.nu,?
            if VERBOSE:
                print("%s Cross Section Saved"%molecule)
        
    @opt.timeit
    def load_HITRAN_single(self,molecule,cache=None):
        """
        Loading molecular cross section for given pressure and temperature
        Will have to be optimized when doing tp profile retrieval.
        This will work when only retrieving on mixing ratio.
        Don't even need to interpolate if matched temperature and pressure
        will see how long the interpolation takes, then judge.
        """
        if cache != None:
            self.xsec[molecule] = cache
    
        xsec = h5py.File("%s/%s.hdf5"%(self.DB_DIR,molecule), "r")
        
        
        self.xsec[molecule] = self.grid_interpolate(np.array(xsec["results"]))
        
    def load_HITRAN_raw_grid(self,molecule):
        
        xsec = h5py.File("%s/%s.hdf5"%(self.DB_DIR,molecule), "r")
        
        return xsec["results"]

    def load_Exomol(self, molecule):
        
        self.xsec = ["Exomol_%s not implemented"%molecule]
            
    def load_NIST(self, molecule):
        
        x1,y1 = load_NIST_spectra(bio_molecule,["wn","T"],True)
        
        Pref = 10000.
        Tref = 300.        
        nref = Pref/(BoltK*Tref)
        lref = 0.05        

        y1 = np.array(y1)+(1-(np.mean(y1)+np.median(y1))/2)
        y1new = []
        for i in y1:
            if i > 1:
                y1new.append(1)
            else:
                y1new.append(i)
        y1 = y1new     
        
        # interpolation
        yinterp = np.interp(self.nu,x1,y1)
        xsec = -np.log(yinterp)/(nref*lref)*10000  # multiply by a factor of 10000 due to unit conversion
    
        
        
        return xsec 

    def load_PNNL(self, molecule):
        
        datapath = "../../SEAS_Input/Cross_Section/HITRAN_Xsec/isoprene/C5-H8_298.1K-760.0K_600.0-6500.0_0.11_N2_505_43.xsc"
        
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
        yinterp = np.interp(self.nu,xlist,ylist)
        
        return yinterp
        
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
        cloud_deck          = float(self.user_input["Xsec"]["Cloud"]["Deck"])
        cloud_opacity       = float(self.user_input["Xsec"]["Cloud"]["Opacity"])
        normalized_pressure = np.array(self.user_input["Prototype"]["Normalized_Pressure"],dtype=float)
        
        #c = Simple_Gray_Cloud_Simulator(cloud_deck,cloud_opacity)
        
        for i,P in enumerate(normalized_pressure):
            if P < cloud_deck:
                normalized_cloud_xsec.append(np.zeros(len(self.nu)))
            else:
                normalized_cloud_xsec.append(np.ones(len(self.nu))*cloud_opacity)
        print("Gray Cloud Loaded")
        return normalized_cloud_xsec   

    def load_mie_cloud(self):
        
        normalized_cloud_xsec = []
        normalized_pressure = np.array(self.user_input["Prototype"]["Normalized_Pressure"],dtype=float)
        
        # will update this to np version once testing is done
        cloud_sigma = self.calculate_mie_cloud()
        for i,P in enumerate(normalized_pressure):
            normalized_cloud_xsec.append(cloud_sigma)
        print("Mie Cloud Loaded")
        
        return normalized_cloud_xsec   

    def calculate_mie_cloud(self):
        
        Source = self.user_input["Xsec"]["Cloud"]["Source"]
        
        lam, n, k = np.genfromtxt("../../SEAS_Input/Refractive_Index/%s"%Source).T
        
        # The miepython module takes the imaginary as negative
        m = n - 1j*k
        
        mean_radius = float(self.user_input["Xsec"]["Cloud"]["Mean_Radius"])
        std         = float(self.user_input["Xsec"]["Cloud"]["Standard_Deviation"])
        sample      = int(self.user_input["Xsec"]["Cloud"]["Sample_Size"])
        
        # draws 10 sample. Not much difference between 10 and 100. 
        # if use 100, may want to output result to avoid duplicated calculations
        radii = np.random.normal(mean_radius,std,sample)
        
        extinct_xsec = np.zeros((len(radii),len(lam)))
        for i,radius in enumerate(radii):
            x = 2*np.pi*radius/lam
            qext, qsca, qback, g = mp.mie(m,x)
            extinct_xsec[i] = qext*np.pi*(radius/1e4)**2 # 1e4 is because um -> cm so that unit is 1/cm^2    
        
        # calculate the cross section given the data
        cloud_sigma = np.mean(np.array(extinct_xsec),axis=0)
        # Return the cross section interpolated to self.nu
        # This only works in linear?
        
        
        return np.interp(10000./self.nu, lam, cloud_sigma)


# Saved for later development
def main_molecule_selector(molecule, 
                           preference_order=["Exomol","HITRAN_Lines","HITRAN_Xsec","NIST"],
                           auxillary=True):
    """
    Due to our data coming from multiple sources, we need a way to automatically   
    select which molecule gets used in the simulation.
    """

    for preference in preference_order:
        
        if preference == "Exomol":
            pass
        elif preference == "HITRAN_Lines":
            pass
        elif preference == "HITRAN_Xsec":
            pass
        elif preference == "NIST":
            pass
            





    
    


