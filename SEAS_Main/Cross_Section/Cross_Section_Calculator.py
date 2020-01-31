"""

Main module for handling lower level operations of cross section generation

This module replaces the old
cross_section_calculator and cross_section_database_generator

Citation for Hapi:
R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016) [Link to article].


"""

import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.External_Utils.hapi as hp

from SEAS_Main.Cross_Section.HITRAN_Match import HITRAN_Match

def hitran_arange_(lower,upper,step):
    npnt = floor((upper-lower)/step)+1
    upper_new = lower + step*(npnt-1)
    if abs((upper-upper_new)-step) < 1e-10:
        upper_new += step
        npnt += 1    
    return linspace(lower,upper_new,npnt)

def calculate_pressure_layers(P_surface = 100000,P_Cutoff = 0.00001):
    """
    Generates pressure as a function of increasing scale height 
    based on surface pressure and pressure cutoff
    
    The number of significant digit is set to 3, but it really doesn't matter
    because the cross section used for spectra simulation 
    will be interpolated from the cross section grid. 
    """
    layers = np.ceil(-np.log(P_Cutoff/P_surface))    
    return [float("%.3g"%x) for x in np.exp(-np.arange(layers))*P_surface]

def calculate_temperature_layers(T_Min=100, T_Max=800, Step=50):
    """
    Generate the temperature layers. +1 so that T_Max is inclusive
    """
    return np.arange(T_Min,T_Max+1,Step)

class cross_section_calculator():

    def __init__(self,d_path,r_path,molecule,component,numin,numax,step=0.1):
        """
        P: Unit is Pa
        n: molecule name
        m: molecular isotope
        i: molecule abundance
        """
        self.d_path = d_path
        self.r_path = r_path
        
        self.molecule = molecule
        
        self.n = component[0]
        self.m = component[1]
        self.i = component[2]

        self.numin = numin
        self.numax = numax   
        self.step  = step

        hp.db_begin(self.d_path)
        hp.fetch(molecule,self.n,self.m,numin,numax)
        #hp.select(molecule,ParameterNames=("nu","sw"), Conditions=("between","nu",float(self.numin),float(self.numax)))
    
    def hapi_calculator(self,P=1.,T=300.,gamma="gamma_self",cross=True):
        """
        This calculation is slow and can not optimize due to external package dependencies, 
        However, this is a "done once" deal so should be ok unless research involve tweeking different cross section types
        This module can be replaced with alternative methods or user developed methods.
        
        issues with the select function in hapi not working as intended
        """
        
        self.P = P
        self.T = T
        self.gamma = gamma
        self.cross = cross
        
        nu, coef = hp.absorptionCoefficient_Voigt(((self.n,self.m,self.i),),
                                                  self.molecule, 
                                                  OmegaStep=self.step,
                                                  HITRAN_units=self.cross,
                                                  GammaL=self.gamma, 
                                                  Environment={'p':float(self.P),'T':float(self.T)})
    
        return nu,coef













        
        
        
        