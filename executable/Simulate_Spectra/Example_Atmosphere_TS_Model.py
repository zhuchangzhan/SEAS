

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
    
    user_input["turnoff"] = ""
    
    simulation = TS.Transmission_Spectra_Simulator(user_input)
    simulation.load_atmosphere_geometry_model()
    
    
    
    return simulation.user_input

@opt.timeit
def Simulate_Atmosphere_Observation(user_input):
       
    nu_               = user_input["Spectra"]["Wavelength"]
    Atmosphere_Height = user_input["Spectra"]["Atmosphere_Height"]
    
    noise = Photon_Noise(user_input)
    bin_edges, bin_width, bin_centers = noise.determine_bin()
    nu,  convolved_Height = noise.calculate_convolve(nu_, Atmosphere_Height)    
    D_atmosphere, bin_edges, binnumber = stats.binned_statistic(10000./nu_[::-1], Atmosphere_Height[::-1], bins=bin_edges)

    user_input["Spectra"]["bin_width"]    = bin_width
    user_input["Spectra"]["bin_centers"]  = bin_centers
    user_input["Spectra"]["nu"]           = nu
    user_input["Spectra"]["bin_values"]   = convolved_Height
    user_input["Spectra"]["D_atmosphere"] = D_atmosphere
    
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
    
    return user_input

@opt.timeit
def cal_binned_SNR(nu, bio_depth,bio_transit_depth_bin,b_SNR_bio,b_ATM_SNR_bio):

    #plt.plot(10000./nu,(bio_depth-ref_depth)*1e6)
    #plt.plot(b_wav,b_SNR_bio)
    #plt.plot(b_wav,b_SNR_ref)
    
    flat_line_a = []
    flat_line_b = []
    flat_line_c = []
    
    for wavs,val in zip(10000./nu,bio_depth):
        if wavs >= 1 and wavs < 5:
            flat_line_a.append(val)
        if wavs >= 5 and wavs < 12:
            flat_line_b.append(val)
        if wavs >= 12:
            flat_line_c.append(val)
    
    mean_atmosphere_radius_per_wavelength = bio_transit_depth_bin
    standard_error_atmosphere_radius_per_wavelength = bio_transit_depth_bin/b_SNR_bio
    
    a,b,c = 33, 18, 15
    
    mean_atmosphere_radius_a = np.average(flat_line_a)
    mean_atmosphere_radius_b = np.average(flat_line_b)
    mean_atmosphere_radius_c = np.average(flat_line_c)
    
    mean_atmosphere_radius_ = np.array([mean_atmosphere_radius_a,
                               mean_atmosphere_radius_b,
                               mean_atmosphere_radius_c])
                                        
    
    standard_error_atmosphere_radius_ = np.array(mean_atmosphere_radius_)/b_ATM_SNR_bio
    
    # extrapolated for subtraction
    mean_atmosphere_radius = np.concatenate([mean_atmosphere_radius_a*np.ones(a),
                                             mean_atmosphere_radius_b*np.ones(b),
                                             mean_atmosphere_radius_c*np.ones(c)])
    ATM_SNR_a,ATM_SNR_b,ATM_SNR_c = b_ATM_SNR_bio
    ATM_SNR = np.concatenate([ATM_SNR_a*np.ones(a),
                              ATM_SNR_b*np.ones(b),
                              ATM_SNR_c*np.ones(c)])
    
    standard_error_atmosphere_radius = mean_atmosphere_radius/ATM_SNR
    
    
    atm_diff = mean_atmosphere_radius_per_wavelength-mean_atmosphere_radius
    atm_error = np.sqrt(standard_error_atmosphere_radius_per_wavelength**2+standard_error_atmosphere_radius**2)
    
    
    Max_SNR = np.max(abs((atm_diff)/atm_error))
    
    
    return Max_SNR,mean_atmosphere_radius_per_wavelength,standard_error_atmosphere_radius_per_wavelength,mean_atmosphere_radius_,standard_error_atmosphere_radius_,atm_diff,atm_error
    
@opt.timeit
def display_output_spectra(user_input,file):
    
    
    fig = plt.figure(figsize=(12, 6)) 
    ax1 = fig.add_subplot(111) 
    
    R_Planet         = float(user_input["Planet"]["R_Planet"])*R_Earth
    R_Star           = float(user_input["Star"]["R_Star"])*R_Sun
    
    
    nu,flux = user_input["Spectra"]["nu"], user_input["Spectra"]["bin_values"]
    
    ref_depth = (flux+R_Planet)**2/R_Star**2

    #ax1.plot(10000./nu,ref_depth*1e6,label = file.split("/")[-1])
        
    ax1.plot(10000./nu,(flux+R_Planet)**2/R_Star**2,label = file.split("/")[-1])
    
    """
    with open("../Pandexo_Integration/%s_%s"%(file.split("/")[1],file.split("/")[-1].replace(".dat",".txt")),"w") as f:
        for i,j in zip(10000./nu[::-1],ref_depth[::-1]):
            f.write("%s %s\n"%(i,j))
    """
    
    #ax1.plot(10000./user_input["Xsec"]["nu"],user_input["Xsec"]["Cloud"]["Value"][0],label = mean)

    ax1.set_ylabel("Transit Depth (ppm)", fontsize=16) 
    #ax1.set_xlabel('Wavelength ($\mu$m)', fontsize=22)
    ax1.legend(fontsize=10)
    ax1.set_xlim(0.6,25)
    ax1.grid(alpha=0.3)
    ax1.set_xscale("log")
    ax1.axvspan(0.6,1,facecolor="g",alpha=0.2)
    ax1.axvspan(1,5,facecolor="g",alpha=0.4)
    ax1.axvspan(5,12,facecolor="r",alpha=0.4)
    ax1.axvspan(12,25,facecolor="r",alpha=0.2)
    #ax1.axhline(observed_radius*1e6,color="0.5")
    ax1.set_xticks([1,2,3,4,5,6,7,8,9,10,12,15,20,25])
    ax1.set_xticklabels(["1","2","3","4","5","6","7","8","9","10","12","15","20","25"])
    ax1.tick_params(which='both',direction="in",labelsize=16,
                labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                bottom=True, top=True, left=True, right=True)  
    ax1.set_xlabel('Wavelength ($\mu$m)', fontsize=16)
    
    plt.show()
    
@opt.timeit
def Forward_Model_Architecture():
    
    global user_input
    user_input["Prototype"]["Source"]="Photochemistry"
    file = os.path.join("../../SEAS_Input/Atmosphere_Data/Atmosphere_Prototype/Example",
                         user_input["Prototype"]["Scenario_File"])
    
    
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
    
    
    display_output_spectra(user_input,file)
    

if __name__ == "__main__":
    
    
    
    Forward_Model_Architecture()
















