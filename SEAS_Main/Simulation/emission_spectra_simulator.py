
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
        
        T_Surface = normalized_temperature[0]
        
        Total_Intensity = np.zeros(len(nu))
        Surface_Intentsity = calc.blackbody_nu(nu,T_Surface) 
        
        tau_of_each_layer = []
        layer_count = np.arange(len(normalized_pressure))
        
        for C,T,P,L in zip(layer_count,normalized_temperature,normalized_pressure,normalized_scale_height):
            
            tau_per_layer = np.zeros(len(nu))
            
            for moi, molecule in enumerate(normalized_molecules):
                sigma_cm = normalized_cross_section[molecule][C]
                sigma = sigma_cm*0.0001
                n = P/(k*T)*normalized_abundance[molecule][C]
                tau_per_layer+=n*sigma*L
                
            Surface_Intentsity *= np.exp(-tau_per_layer)   
            tau_of_each_layer.append(tau_per_layer)
            
        Total_Intensity+=Surface_Intentsity
            
        for C,T,P,L in zip(layer_count,normalized_temperature,normalized_pressure,normalized_scale_height):
            
            tau_cur_layer = tau_of_each_layer[C]
            Layer_Intentsity = calc.blackbody_nu(nu,T)*tau_cur_layer
            for j in layer_count[C:]:
                Layer_Intentsity *= np.exp(-tau_of_each_layer[j])
            Total_Intensity += Layer_Intentsity    
        
        
        self.user_input["Spectra"]["Wavelength"]      = nu
        self.user_input["Spectra"]["Total_Intensity"] = Total_Intensity
        
    