
import os
import sys
import numpy as np
    
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Main.Physics.astrophysics as calc
import matplotlib.pyplot as plt

class Transmission_Spectra_Simulator():
    
    def __init__(self,user_input):
        self.user_input = user_input

    def load_boxcar_model(self):
        
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = float(prototype["Normalized_Pressure"][0])
        normalized_temperature      = float(prototype["Normalized_Temperature"][0])
        normalized_scale_height     = float(prototype["Normalized_Scale_Height"][0])
        boxcar_pathlength           = 1
        
        normalized_molecules        = prototype["Molecule_List"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_rayleigh         = self.user_input["Xsec"]["Rayleigh"]["Value"]

        ChunkTau = np.zeros(len(nu))      
        
        
        for molecule in normalized_molecules:  
            molecular_ratio = normalized_abundance[molecule][0]
            number_density  = (normalized_pressure/(BoltK*normalized_temperature))*molecular_ratio
            sigma           = normalized_cross_section[molecule][0]

            # 0.0001 is convertion factor for xsec from cgs to SI
            ChunkTau += number_density*sigma*boxcar_pathlength*0.0001 
        
        BeamTrans = calc.calc_transmittance(ChunkTau) 
        Atmosphere_Height = (1-BeamTrans)*normalized_scale_height
            
        
        self.user_input["Spectra"]["Wavenumber"] = nu
        self.user_input["Spectra"]["Total_Transit_Signal"] = BeamTrans
        
    def load_boxcar_model_multilayer(self):
        
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"])
        normalized_temperature      = np.array(prototype["Normalized_Temperature"])
        normalized_scale_height     = np.array(prototype["Normalized_Scale_Height"])
        
        boxcar_pathlength           = [1,1,1]
        
        normalized_molecules        = prototype["Molecule_List"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_rayleigh         = self.user_input["Xsec"]["Rayleigh"]["Value"]

       
        layers = len(normalized_pressure)
        
        Atmosphere_Height = np.zeros(len(nu))      
        
        for cur in range(layers): # assuming each beam has only 1 layer
            
            BeamTau = np.zeros(len(nu))      
            
            p_cur = normalized_pressure[cur]
            t_cur = normalized_temperature[cur]
            l_cur = boxcar_pathlength[cur]
            
            ChunkTau = np.zeros(len(nu))      
             
            for molecule in normalized_molecules:  
                molecular_ratio = normalized_abundance[molecule][cur]
                number_density  = (p_cur/(BoltK*t_cur))*molecular_ratio*1e-6 # molc/m^3 -> molc/cm^3
                sigma           = normalized_cross_section[molecule][cur] #cm^2/molc
    
                # 0.0001 is convertion factor for xsec from cgs to SI
                # cm^2/molc * molc/cm^3 -> cm^-1
                ChunkTau += number_density*sigma*l_cur*100
                
                
            # Sum up all the absorption along the beam 
            BeamTau += ChunkTau   
            
            BeamTrans = calc.calc_transmittance(BeamTau) 
            effective_height = (1-BeamTrans)*normalized_scale_height[cur]
            
            Atmosphere_Height += effective_height
        
        self.user_input["Spectra"]["Wavenumber"] = nu
        self.user_input["Spectra"]["Atmosphere_Height"]    = Atmosphere_Height
        self.user_input["Spectra"]["Total_Transit_Signal"] = Atmosphere_Height # proxy this for now for retrieval implementation, will change
                              
    def load_atmosphere_geometry_model(self):
        
        # Loading all the respective data from dictionary into parameters
        prototype = self.user_input["Prototype"]
        
        Molecule  = self.user_input["Xsec"]["Molecule"]["Enable"] == "True"
        Cloud     = self.user_input["Xsec"]["Cloud"]["Enable"] == "True"  # bool("False") returns True so... 
        CIA       = self.user_input["Xsec"]["CIA"]["Enable"]   == "True"
        
        # The parameter for biomolecules need to be better defined when there is no biosignature. 
        # Ie if one just want a atmosphere and not a comparison.
        Bio_Molecules               = prototype["Bio_Molecule_List"]
        Particle_Ratio              = float(self.user_input["Xsec"]["Cloud"]["Particle_Ratio"])
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"],dtype=float)
        normalized_temperature      = np.array(prototype["Normalized_Temperature"],dtype=float)
        normalized_molecules        = prototype["Molecule_List"]  # defined in load_Atmosphere_Profile
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        normalized_scale_height     = prototype["Normalized_Scale_Height"] 
        normalized_mean_mw          = prototype["normalized_mean_mw"]    
        
        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_CIA              = self.user_input["Xsec"]["CIA"]
        normalized_rayleigh         = self.user_input["Xsec"]["Rayleigh"]["Value"]
        normalized_cloud_xsec       = self.user_input["Xsec"]["Cloud"]["Value"]
        CIA_species                 = ["H2-H2"]
        
        
        """
        for molecule in normalized_molecules:
            plt.plot(10000./nu,normalized_cross_section[molecule][0])
            plt.yscale("log")
            plt.xscale("log")
        plt.show()
        """

        
        TotalBeams           = len(normalized_pressure)
        Total_Tau            = np.zeros(len(nu))
        Offset               = 0#float(self.user_input["Atmosphere_Effects"]["Base_Line"]["offset"])
        Total_Transit_Signal = np.ones(len(nu))*(self.user_input["Atmosphere"]["Base_TS_Value"]+Offset)
        Total_Height         = np.ones(len(nu))*float(self.user_input["Planet"]["R_Planet"])*R_Earth
        Atmosphere_Height    = np.zeros(len(nu))
        base_layer           = float(self.user_input["Planet"]["R_Planet"])*R_Earth


        # this need to be better implemented
        if self.user_input["Config"]["molecule_turnoff"] == None:
            pass
        else:
            normalized_cross_section[self.user_input["Config"]["molecule_turnoff"]] = np.zeros((TotalBeams,12000))
        

        
        # this need to be parameterized
        Cloud = False
        for i in range(TotalBeams):
            
            if i == 0:
                prev_layer = base_layer
                base_layer += normalized_scale_height[i]
                # skip the bottom layer? this is the beam that "touch" the surface
                # need to think about this when rounding
                continue
    
            # opacity per beam
            BeamTau = []
            prev_pathl = 0
            target_layer = base_layer    
            
            for j in range(TotalBeams-i):
                cur = j+i
                
                target_layer += normalized_scale_height[cur]
                pathl = np.sin(np.arccos(base_layer/target_layer))*target_layer - prev_pathl # m
                prev_pathl += pathl   
                
                mean_mw = normalized_mean_mw[cur]
                P_cur = normalized_pressure[cur]
                T_cur = normalized_temperature[cur]  
                PKT   = P_cur/(BoltK*T_cur)
                
                # opacity per chunk of the beam, this can be thought as the test tube case
                # In the future this part will be spinned out as its seperate function
                # So that it's shared between each type of models
                ChunkTau = []   
                haze_source_mr = {}     
                if Molecule:
                    for molecule in normalized_molecules:    
                        if molecule in Bio_Molecules and Cloud:
                            haze_source_mr[molecule] = normalized_abundance[molecule][cur]*Particle_Ratio
                        
                        #weird how abundance and cross section are wired differently
                        molecular_ratio = normalized_abundance[molecule][cur]
                        number_density = PKT*molecular_ratio*1e-6 # molc/m^3 -> molc/cm^3
                        rayleigh = normalized_rayleigh[molecule]*molecular_ratio
                        sigma = normalized_cross_section[molecule][cur] #cm^2/molc
                        
                        # 2 because mirror other half of the atmosphere
                        # 0.0001 is convertion factor for xsec from cgs to SI
                        # cm^2/molc * molc/cm^3 -> cm^-1
                        ChunkTau_Per_Molecule = number_density*(sigma+rayleigh)*pathl*100*2
            
                        if ChunkTau == []:
                            ChunkTau = ChunkTau_Per_Molecule
                        else:
                            ChunkTau += ChunkTau_Per_Molecule 
                    
                
                if CIA:
                    for mol1mol2 in CIA_species:
                        molecule1, molecule2 = mol1mol2.split("-")
                        
                        molecular_ratio1 = normalized_abundance[molecule1][cur]
                        molecular_ratio2 = normalized_abundance[molecule2][cur]
                        
                        number_density1 = PKT*molecular_ratio1*1e-6
                        number_density2 = PKT*molecular_ratio2*1e-6
                        
                        CIA_sigma = normalized_CIA[mol1mol2][cur] #cm^2/molc
                        
                        Chunk_CIA_Tau = number_density1*number_density2*CIA_sigma*pathl*100*2
                        
                        ChunkTau += Chunk_CIA_Tau      
                    
                # Simulating cloud
                if Cloud:
                    
                    for molecule in Bio_Molecules:
                    
                        # Load cloud particle number density
                        air_rho     = calc.calculate_air_density(P_cur,T_cur,mean_mw)*1e-6  #m^3 to cm^3
                        haze_ratio  = haze_source_mr[molecule] # this is a unitless scaling factor
                        haze_rho    = float(self.user_input["Xsec"]["Cloud"]["Particle_Density"]) # g/cm^3
                        haze_radius = float(self.user_input["Xsec"]["Cloud"]["Mean_Radius"])*1e-4 # radi from um to cm
                        particle_number_density = calc.calc_cloud_number_density(air_rho,haze_ratio,haze_rho,haze_radius)  # particle/cm^3
                        
                        
                        # load cloud cross section
                        cloud_sigma = normalized_cloud_xsec[cur] # cm^2/particle
                        
                        # cm^2/particle * particle/cm^3 -> cm^-1
                        ChunkTau_Per_Cloud = cloud_sigma*particle_number_density*pathl*100*2
                        
                        ChunkTau += ChunkTau_Per_Cloud
                    
                    # grey cloud to be added back later
                    #ChunkTau+= normalized_cloud_xsec[cur]
                    
                # Sum up all the absorption along the beam 
                if BeamTau == []:
                    BeamTau = ChunkTau
                else:
                    BeamTau += ChunkTau   
            
            # This part need to be redone
            BeamTrans = calc.calc_transmittance(BeamTau) 
            effective_height = (1-BeamTrans)*normalized_scale_height[i]
            RingArea = (base_layer**2-prev_layer**2)/(float(self.user_input["Star"]["R_Star"])*R_Sun)**2
            Ring_Transit_Signal = (1-BeamTrans)*RingArea
            
            Total_Transit_Signal += Ring_Transit_Signal
            Atmosphere_Height += effective_height
            
            # update to the next beam up
            prev_layer = base_layer
            base_layer += normalized_scale_height[i]        
        
        self.user_input["Spectra"]["Wavelength"]           = nu
        self.user_input["Spectra"]["Atmosphere_Height"]    = Atmosphere_Height
        self.user_input["Spectra"]["Total_Transit_Signal"] = Total_Transit_Signal
        
        
"""
# fix or move this to somewhere else
for molecule in normalized_molecules:
    if molecule != "HNO3":
        normalized_cross_section[molecule] = np.zeros((24,12000))

for molecule in normalized_rayleigh:
    normalized_rayleigh[molecule] = np.zeros(12000)
    
import matplotlib.pyplot as plt
for i in normalized_cross_section["HNO3"]:
    plt.plot(i)
plt.yscale("log")
plt.show()
"""        
        
        
        
        