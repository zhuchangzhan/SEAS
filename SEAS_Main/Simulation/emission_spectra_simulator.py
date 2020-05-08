
import os
import sys
import numpy as np
    
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Main.Physics.astrophysics as calc

import matplotlib.pyplot as plt

class Emission_Spectra_Simulator():
    
    def __init__(self,user_input):
        self.user_input = user_input
    
    def load_atmosphere_geometry_model(self):
        
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"],dtype=float)
        normalized_temperature      = np.array(prototype["Normalized_Temperature"],dtype=float)
        normalized_scale_height     = prototype["Normalized_Scale_Height"] 
        normalized_molecules        = prototype["Molecule_List"]

        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        
        turnoff = [self.user_input["Config"]["molecule_turnoff"]]
        
        def calc_tau(C,T,P,L):
            
            tau_cur_layer = np.zeros(len(nu))
            
            for moi, molecule in enumerate(normalized_molecules):
                if molecule in turnoff:
                    tau_cur_layer+= np.zeros(len(nu))
                else:
                    sigma_cm = normalized_cross_section[molecule][C]
                    sigma = sigma_cm*0.0001
                    n = P/(k*T)*normalized_abundance[molecule][C]
                    
                    tau_cur_layer+=n*sigma*L
                    
            return tau_cur_layer
            
        # create a 25 x 12000 array of zeros
        Layer_Intensity = np.zeros((len(normalized_pressure)+1,len(nu)))     
        
        # calculate emission (spectral radiance) from surface
        Layer_Intensity[0] = calc.blackbody_nu(nu,normalized_temperature[0]) 
        
        # for each layer of the atmosphere
        for C,(T,P,L) in enumerate(zip(normalized_temperature,
                                       normalized_pressure,
                                       normalized_scale_height)): 
            
            # calculate opacity of the current layer
            tau_cur_layer = calc_tau(C,T,P,L)
            
            # calculate emission (spectral radiance) for the current layer
            Layer_Intensity[C+1] = calc.blackbody_nu(nu,normalized_temperature[C])*(1-np.exp(-tau_cur_layer))
            
            # attenuation of each layer below by current layer 
            for j in range(C+1):
                Layer_Intensity[j] *= np.exp(-tau_cur_layer)        
        
            
        Total_Intensity = np.sum(Layer_Intensity,axis=0) 
        
        self.user_input["Spectra"]["Wavelength"]      = nu
        self.user_input["Spectra"]["Total_Intensity"] = Total_Intensity      
    
    
        
    
    