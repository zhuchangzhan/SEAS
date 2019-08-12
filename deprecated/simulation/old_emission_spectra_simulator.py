#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Emission Spectra Simulator
"""


import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Main.simulation.transmission_spectra_simulator import TS_Simulator
from scipy.constants import h,k,c
from SEAS_Utils.common_utils.constants import *

class ES_Simulator(TS_Simulator):
    
    
    def __init__(self, user_input):
        
        TS_Simulator.__init__(self, user_input)
    
    def blackbody_lam(self, lam, T):
        """ Blackbody as a function of wavelength (um) and temperature (K).
    
        returns units of erg/s/cm^2/cm/Steradian
        """
    
        lam = 1e-6 * lam # convert to metres
        return 2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))


    def load_atmosphere_geometry_model(self):
        

        
        Total_Layers = len(self.normalized_pressure)
        
        normalized_pressure         = self.normalized_pressure
        normalized_temperature      = self.normalized_temperature
        normalized_molecules        = self.normalized_molecules
        normalized_abundance        = self.normalized_abundance
        normalized_cross_section    = self.normalized_cross_section
        normalized_scale_height     = self.normalized_scale_height 
        
        

        Surface_T = normalized_temperature[0]
        Base_EM = self.blackbody_lam(10000./self.nu, Surface_T)
        
        
        Total_EM = Base_EM
        

        for i in range(Total_Layers):
        
            for m, molecule in enumerate(normalized_molecules):
                
                pathl = normalized_scale_height[i]
                sigma = normalized_cross_section[m][i][i]
                molecular_ratio = normalized_abundance[i][m]
                n = (normalized_pressure[i]/(BoltK*normalized_temperature[i]))*molecular_ratio    
                        
                tau = n*sigma*pathl*0.0001
                
        
                Layer_Transmittance = np.e**(-tau)

                Layer_Absorbance = 1-Layer_Transmittance 
            
            
                Layer_TS = Total_EM*Layer_Transmittance
                
                Layer_EM = self.blackbody_lam(10000./self.nu, normalized_temperature[i])*Layer_Absorbance
                
                Total_EM = Layer_TS + Layer_EM
            
        
        
        
        
        
        
        
        return Total_EM#*np.e**(-tau)










