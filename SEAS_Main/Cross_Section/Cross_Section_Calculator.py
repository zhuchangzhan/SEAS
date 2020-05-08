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
    Helper function for generating temperature layers
    +1 is added to T_Max so that it is inclusive [100,800]
    """
    return np.arange(T_Min,T_Max+1,Step)

class cross_section_calculator():

    def __init__(self,d_path,molecule,component,numin,numax,step=0.1,remake=False):
        """
        Initiate cross section calculator instance with required parameters
        
        Parameters
        ----------
        d_path : str
            filepath for where the linelist data will be stored
        molecule : str
            Name of the molecule. Has to be one which HITRAN recognize
        numin : float
            minimum wavenumber [cm^-1]
        numax : float
            maximum wavenumber [cm^-1]
        step : float
            wavenumber increment step size. 
        remake : bool
            Flag for whether the hitran linelist will be downloaded again or not. 
            Duplicate with SL_Flag. Will remove in future iterations but keep now for compat.

        """
        self.d_path = d_path
        
        self.molecule = molecule
        
        self.n = component[0]
        self.m = component[1]
        self.i = component[2]

        self.numin = numin
        self.numax = numax   
        self.step  = step

        hp.db_begin(self.d_path)
        
        
        if os.path.isfile(os.path.join(self.d_path,"%s.header"%(molecule))) and not remake:
            pass
            #hp.select(molecule)#,ParameterNames=("nu","sw"), Conditions=("between","nu",float(self.numin),float(self.numax)))
        else:
            print("getting new data")
            hp.fetch(molecule,self.n,self.m,numin,numax)

            
    def hapi_calculator(self,P=1.,T=300.,gamma="gamma_self",cross=True):
        """
        Calculate cross section using the hapi package.
        
        This calculation is slow and difficult to optimize due to external package dependencies, 
        However, cross section grid generation should be an "onetime" task for most user applications.
        
        Parameters
        ----------        

        P : float
            Pressure [bar]
        T : float
            Temperature [K]        
        gamma : str
            Hapi parameter 
        cross : bool
            If True, use cm^2/molecule for output cross section value
            This should be kept as True unless you know what you're doing
        
        Returns
        -------
        nu : array
            wavenumber array
        coef : array
            cross section array
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
                                                  Environment={'p':P,'T':T})
    
        return nu,coef













        
        
        
        