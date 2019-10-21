import os
import sys
import numpy as np
    
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Main.Physics.astrophysics as calc


class Transmission_Spectra_Simulator():
    
    def __init__(self,user_input):
        self.user_input = user_input

    def load_boxcar_room_model(self):
        
        prototype = self.user_input["Prototype"]
        
        boxcar_pressure             = 100000
        boxcar_temperature          = 300
        boxcar_pathlength           = 1
        
        normalized_molecules        = prototype["Molecule_List"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_rayleigh         = self.user_input["Xsec"]["Rayleigh"]["Value"]

        ChunkTau = np.zeros(len(nu))        
        for molecule in normalized_molecules:       
            
            molecular_ratio = normalized_abundance[molecule]
            number_density  = (normalized_pressure/(BoltK*normalized_temperature))*molecular_ratio
            sigma           = normalized_cross_section[molecule]

            # 0.0001 is convertion factor for xsec from cgs to SI
            ChunkTau += number_density*sigma*boxcar_pathlength*0.0001 

        return nu,calc.calc_transmittance(ChunkTau) 
        
    def load_atmosphere_geometry_model(self):
        
        # Loading all the respective data from dictionary into parameters
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"],dtype=float)
        normalized_temperature      = np.array(prototype["Normalized_Temperature"],dtype=float)
        normalized_molecules        = prototype["Molecule_List"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        normalized_scale_height     = prototype["Normalized_Scale_Height"]     
        
        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_rayleigh         = self.user_input["Xsec"]["Rayleigh"]["Value"]
        normalized_cloud_xsec       = self.user_input["Xsec"]["Cloud"]["Value"]
        
        TotalBeams           = len(normalized_pressure)
        Total_Tau            = np.zeros(len(nu))
        Offset               = 0#float(self.user_input["Atmosphere_Effects"]["Base_Line"]["offset"])
        Total_Transit_Signal = np.ones(len(nu))*(self.user_input["Atmosphere"]["Base_TS_Value"]+Offset)
        Total_Height         = np.ones(len(nu))*float(self.user_input["Planet"]["R_Planet"])*R_Earth
        Atmosphere_Height    = np.zeros(len(nu))
        base_layer           = float(self.user_input["Planet"]["R_Planet"])*R_Earth
        
        Cloud, CIA = True,False
        
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
                pathl = np.sin(np.arccos(base_layer/target_layer))*target_layer - prev_pathl
                prev_pathl += pathl   
                
                
                # opacity per chunk of the beam, this can be thought as the test tube case
                ChunkTau = []        
                for molecule in normalized_molecules:       
                    
                    #weird how abundance and cross section are wired differently
                    molecular_ratio = normalized_abundance[molecule][cur]
                    number_density = (normalized_pressure[cur]/(BoltK*normalized_temperature[cur]))*molecular_ratio
                    rayleigh = normalized_rayleigh[molecule]*molecular_ratio
                    sigma = normalized_cross_section[molecule][cur]
        
                    # 2 because mirror other half of the atmosphere
                    # 0.0001 is convertion factor for xsec from cgs to SI
                    ChunkTau_Per_Molecule = number_density*(sigma+rayleigh)*pathl*2*0.0001
        
                    if ChunkTau == []:
                        ChunkTau = ChunkTau_Per_Molecule
                    else:
                        ChunkTau += ChunkTau_Per_Molecule 
                
                
                

                """
                if CIA:
                    for k,CIA_data in enumerate(normalized_CIA):
                        
                        molecule1, molecule2 = self.CIA_File[k].replace("_","-").split("-")[:2]
                        
                        molecular_ratio1 = normalized_abundance_ref[molecule1][j+i]
                        molecular_ratio2 = normalized_abundance_ref[molecule2][j+i]
                        
                        number_density1 = (normalized_pressure[j+i]/(BoltK*normalized_temperature[j+i]))*molecular_ratio1
                        number_density2 = (normalized_pressure[j+i]/(BoltK*normalized_temperature[j+i]))*molecular_ratio2
                        
                        CIA_sigma = CIA_data[j+i]
                        Chunk_CIA_Tau = number_density1*number_density2*CIA_sigma*pathl*2*100*100**-3*100**-3
                        ChunkTau += Chunk_CIA_Tau      
                """
                # Simulating cloud
                if Cloud:
                   ChunkTau+= normalized_cloud_xsec[cur]
                    
                # Sum up all the absorption along the beam 
                if BeamTau == []:
                    BeamTau = ChunkTau
                else:
                    BeamTau += ChunkTau   
            
            BeamTrans = calc.calc_transmittance(BeamTau) 
            effective_height = (1-BeamTrans)*normalized_scale_height[i]
            RingArea = (base_layer**2-prev_layer**2)/(float(self.user_input["Star"]["R_Star"])*R_Sun)**2
            Ring_Transit_Signal = (1-BeamTrans)*RingArea
            
            Total_Transit_Signal += Ring_Transit_Signal
            Atmosphere_Height += effective_height
            
            
            # update to the next beam up
            prev_layer = base_layer
            base_layer += normalized_scale_height[i]        
        
        self.user_input["Spectra"]["Wavelength"] = nu
        self.user_input["Spectra"]["Atmosphere_Height"] = Atmosphere_Height
        self.user_input["Spectra"]["Total_Transit_Signal"] = Total_Transit_Signal
        
        
        