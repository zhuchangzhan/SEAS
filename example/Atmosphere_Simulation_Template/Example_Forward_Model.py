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


def blackbody_lam(wav, T):
    """ Blackbody as a function of wavelength (m) and temperature (K).
    returns units of erg/s/cm^2/cm/Steradian
    """
    a = 2*HPlanck*CLight**2
    b = HPlanck*CLight/(wav*BoltK*T)
    intensity = a/((wav**5)*(np.exp(b)-1.0))
    return intensity

@opt.timeit
def Generate_Atmosphere_Spectra(user_input):
    
    simulation = TS.Transmission_Spectra_Simulator(user_input)
    simulation.load_atmosphere_geometry_model()
    
    return simulation.user_input

@opt.timeit
def Simulate_Atmosphere_Observation(user_input):
    
    def sims(user_input):
        nu_            = user_input["Spectra"]["Wavelength"]
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
        
    
        return user_input,nu,convolved_Height,bin_centers,bin_width,D_atmosphere

    def simulated_observation(molecule,stellar,nu,bio_trans,bin_centers,bin_width,D_atmosphere,observe_hour):
    
        if stellar == "Sun":
            T_Star = 5770.
            R_Star = R_Sun
            
        elif stellar == "Ma":
            T_Star = 3000.
            R_Star = 0.26*R_Sun
        
        elif stellar == "Mq":
            T_Star = 3000.
            R_Star = 0.26*R_Sun
                
        #bio_molecule = "PH3"
        #molecule = "CO2"
        #stellar = "Sun"
        g_atm = 3200. 
        if molecule == "H2":
            mu_atm = 4.6*amu2g
        elif molecule == "N2":
            mu_atm = 28*amu2g
        elif molecule == "CO2":
            mu_atm = 44*amu2g   
    
        R_obs        = 6.5
        R_Planet     = 6400000*1.75
        Distance     = 10*Psec  # 10 pc to meter
        Duration     = observe_hour*3600. # 10 hr to second
        Quantum      = 0.3 
        Noise_M      = 1.5
        
        # calculate number of photons 
        B_Body = blackbody_lam(bin_centers*10**-6, T_Star)
        Bin_width = bin_width*10**-6
        A_Star = np.pi*R_Star**2
        Psi_Tele = np.pi*R_obs**2/Distance**2
        E_Total = B_Body*Bin_width*A_Star*Psi_Tele*Duration
        star_photon = (E_Total*bin_centers*10**-6)/(HPlanck*CLight)*Quantum
        
        #signal = (2*R_planet*D_atmosphere)/R_Star**2*num_photon
        
        R_Planet = R_Planet+D_atmosphere
        
        signal = R_Planet**2/R_Star**2*star_photon
        photon_noise = Noise_M*np.sqrt(star_photon)
        
        
        SNR = signal/photon_noise
        
        a_sig,b_sig,c_sig = [],[],[]
        a_noise,b_noise,c_noise = [],[],[]
        for i,sig,noise in zip(bin_centers,signal,photon_noise):
            if i> 1 and i < 5:
                a_sig.append(sig)
                a_noise.append(noise)
            if i >= 5 and i < 12:
                b_sig.append(sig)
                b_noise.append(noise)
            if i >= 12:
                c_sig.append(sig)
                c_noise.append(noise)
        SNR_ATM = [np.sum(a_sig)/np.sqrt(np.sum(np.array(a_noise)**2)),
                   np.sum(b_sig)/np.sqrt(np.sum(np.array(b_noise)**2)),
                   np.sum(c_sig)/np.sqrt(np.sum(np.array(c_noise)**2))]
        
        #SNR_ATM = np.sum(signal)/np.sqrt(np.sum(photon_noise**2))
        
        
        depth = R_Planet**2/R_Star**2#*1e6
        
        
        #bin_means_bio,SNR_bio = simulate_observation(noise,bin_edges,nu_,Bio_Transit_Signal)
        
        return R_Planet,depth,SNR,SNR_ATM
        
    def cal_SNR(bio_depth,bio_transit_depth_bin,b_SNR_bio,b_ATM_SNR_bio):

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

    user_input,nu,ref_trans,b_wav,b_width,b_D_atm_ref = sims(user_input)
    
    main_mol = "CO2"
    stellar = "Ma"
    observe_hour = 100
    R_Planet = 6400000
    R_Star = 0.26*R_Sun
    
    b_R_planet_ref,b_depth_ref,b_SNR_ref,b_ATM_SNR_ref = simulated_observation(main_mol,stellar,nu,ref_trans,b_wav,b_width,b_D_atm_ref,observe_hour)

    ref_depth = (ref_trans+R_Planet)**2/R_Star**2
    ref_depth_bin = (b_D_atm_ref+R_Planet)**2/R_Star**2

    result = cal_SNR(ref_depth,ref_depth_bin,b_SNR_ref,b_ATM_SNR_ref)
    ATM_Max_SNR_ref,aaa,bbb,ccc,ddd,atm_diff_ref,atm_error_ref = result

    mean_atmosphere_radius_per_wavelength_ref = aaa
    standard_error_atmosphere_radius_per_wavelength_ref = bbb
    mean_atmosphere_radius_ref_ = ccc
    standard_error_atmosphere_radius_ref_ = ddd 

    mean_atmosphere_wav     = np.array([3,8.5,18.5])
    mean_atmosphere_width   = np.array([2,3.5,6.5])

    fig = plt.figure(figsize=(12, 9)) 
    ax1 = fig.add_subplot(111) 
    ax1.plot(10000./nu,ref_depth*1e6,label="w/o PH$_3$")



    ax1.errorbar(b_wav,
                 (mean_atmosphere_radius_per_wavelength_ref)*1e6,
                 yerr=standard_error_atmosphere_radius_per_wavelength_ref*1e6,
                 capsize=4, markersize=6,fmt='--o',label="w/o PH$_3$ Binned",color="0.5")
    ax1.errorbar(mean_atmosphere_wav,
                 (mean_atmosphere_radius_ref_)*1e6,
                 xerr=mean_atmosphere_width,
                 yerr=standard_error_atmosphere_radius_ref_*1e6,
                 capsize=4, markersize=6,fmt='o',label="w/o PH$_3$ Atm. Avg.",color="y")   
    

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
def display_output_spectra(user_input):
    
    x = user_input["Spectra"]["Wavelength"]
    y = user_input["Spectra"]["Atmosphere_Height"]
    
    plt.xscale("log")
    plt.plot(10000./x,y)
    plt.show()
    
@opt.timeit
def Forward_Model_Architecture():
    
    global user_input
    
    file = os.path.join("../../SEAS_Input/Atmosphere_Data/Atmosphere_Prototype/Example",
                         user_input["Prototype"]["Scenario_File"])
    
    # Loading TP_Profile, MR_Profile, Molecule List, 
    user_input = load.load_Atmosphere_Profile(user_input, source="Photochemistry", scenario_file=file)
    
    # Hash created for identifying the xsec used for the TP profile
    user_input["Data_IO"]["Hash"] = hashlib.sha224(str(user_input["Prototype"]["TP_Profile"]).encode()).hexdigest()[:8]
    
    # Load Surface_Gravity, Base_TS_Value
    user_input = load.load_Astrophysical_Properties(user_input)
    
    # Load absorption cross section for all molecule and effect (CIA, Cloud, etc)
    user_input = load.load_Absorption_Cross_Section(user_input)
    
    # Load Atmosphere model and generate theoretical spectra
    user_input = Generate_Atmosphere_Spectra(user_input)
    
    # Load Observation parameters and generate simulated observation
    user_input = Simulate_Atmosphere_Observation(user_input)
    
    # Display the output result
    #display_output_spectra(user_input)
    


if __name__ == "__main__":
    
    
    
    Forward_Model_Architecture()
















