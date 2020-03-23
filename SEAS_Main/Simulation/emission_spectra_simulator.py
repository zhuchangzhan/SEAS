
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
    
    def load_atmosphere_geometry_model_old(self):
        
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"],dtype=float)
        normalized_temperature      = np.array(prototype["Normalized_Temperature"],dtype=float)
        normalized_scale_height     = prototype["Normalized_Scale_Height"] 
        normalized_molecules        = prototype["Molecule_List"]

        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        T_Surface = normalized_temperature[0]

        Surface_Intensity = calc.blackbody_nu(nu,T_Surface)
        
        tau_of_each_layer = []
        layer_count = np.arange(len(normalized_pressure))
        
        
        turnoff = self.user_input["turnoff"]
        Total_Intensity = np.zeros(len(nu))
        
        for C,T,P,L in zip(layer_count,normalized_temperature,normalized_pressure,normalized_scale_height):
            
            tau_cur_layer = np.zeros(len(nu))
            
            for moi, molecule in enumerate(normalized_molecules):
                if molecule in turnoff:
                    tau_cur_layer+= np.zeros(len(nu))
                else:
                    sigma_cm = normalized_cross_section[molecule][C]
                    sigma = sigma_cm*0.0001
                    n = P/(k*T)*normalized_abundance[molecule][C]
                    tau_cur_layer+=n*sigma*L
                    
            Surface_Intensity *= np.exp(-tau_cur_layer) 
            tau_of_each_layer.append(tau_cur_layer)
            
        Total_Intensity+=Surface_Intensity
            
        for C,T,P,L in zip(layer_count,normalized_temperature,normalized_pressure,normalized_scale_height):
            
            tau_cur_layer = tau_of_each_layer[C]
            Layer_Intensity = calc.blackbody_nu(nu,T)*tau_cur_layer
            for j in layer_count[C:]:
                Layer_Intensity *= np.exp(-tau_of_each_layer[j])
            Total_Intensity += Layer_Intensity    

        
        
        
        self.user_input["Spectra"]["Wavelength"]      = nu
        self.user_input["Spectra"]["Total_Intensity"] = Total_Intensity

    def load_atmosphere_geometry_model_redo(self):
        
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"],dtype=float)
        normalized_temperature      = np.array(prototype["Normalized_Temperature"],dtype=float)
        normalized_scale_height     = prototype["Normalized_Scale_Height"] 
        normalized_molecules        = prototype["Molecule_List"]

        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        
        turnoff = self.user_input["turnoff"]
        
        
        
        def calc_tau(C,T,P):
            
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
            
        
        layer_count = np.arange(len(normalized_pressure))
        Layer_Intensity = np.zeros((len(normalized_pressure),len(nu)))
        
        Layer_Intensity[0] = calc.blackbody_nu(nu,normalized_temperature[0])

        for C,T,P,L in zip(layer_count,normalized_temperature,normalized_pressure,normalized_scale_height):
            
            tau_cur_layer = calc_tau(C,T,P)
            
            if C == 0:
                Layer_Intensity[C] += calc.blackbody_nu(nu,normalized_temperature[C])*tau_cur_layer
            else:
                Layer_Intensity[C] = calc.blackbody_nu(nu,normalized_temperature[C])*tau_cur_layer
            
            for j in layer_count[:C+1]:
                Layer_Intensity[j] *= np.exp(-tau_cur_layer)
            
        Total_Intensity = np.sum(Layer_Intensity,axis=0)
        
        self.user_input["Spectra"]["Wavelength"]      = nu
        self.user_input["Spectra"]["Total_Intensity"] = Total_Intensity      
    
    def load_atmosphere_geometry_model(self):
        
        prototype = self.user_input["Prototype"]
        
        normalized_pressure         = np.array(prototype["Normalized_Pressure"],dtype=float)
        normalized_temperature      = np.array(prototype["Normalized_Temperature"],dtype=float)
        normalized_scale_height     = prototype["Normalized_Scale_Height"] 
        normalized_molecules        = prototype["Molecule_List"]

        nu                          = self.user_input["Xsec"]["nu"]
        normalized_cross_section    = self.user_input["Xsec"]["Molecule"]
        normalized_abundance        = prototype["Normalized_MR_Profile"]
        
        
        turnoff = self.user_input["turnoff"]
        
        def calc_tau(C,T,P):
            
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
            
        
        layer_count = np.arange(len(normalized_pressure))
        Layer_Intensity = np.zeros((len(normalized_pressure),len(nu)))

        for C,T,P,L in zip(layer_count,normalized_temperature,normalized_pressure,normalized_scale_height):
            
            tau_cur_layer = calc_tau(C,T,P)
            
            if C == 0: # surface
                Layer_Intensity[C] = calc.blackbody_nu(nu,normalized_temperature[C])
            else:
                Layer_Intensity[C] = calc.blackbody_nu(nu,normalized_temperature[C])*(1-np.exp(-tau_prev_layer))
        
            for j in layer_count[:C+1]:
                Layer_Intensity[j] *= np.exp(-tau_cur_layer)
            
            tau_prev_layer = tau_cur_layer
            
            
        Total_Intensity = np.sum(Layer_Intensity,axis=0)
        
        self.user_input["Spectra"]["Wavelength"]      = nu
        self.user_input["Spectra"]["Total_Intensity"] = Total_Intensity      
    
    
        
    
    