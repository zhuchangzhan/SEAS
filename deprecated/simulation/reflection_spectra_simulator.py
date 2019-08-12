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

Reflection Spectra Simulator


Things to consider

excited chromophore => fluorescence
Fresnel equation
refractive index
extinction coefficient

surface and atmosphere reflection and scattering


"""
import os
import sys
import numpy as np
import time
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)

import SEAS_Main.atmosphere_geometry
import SEAS_Main.atmosphere_property
from SEAS_Main.atmosphere_effects.biosig_molecule import load_NIST_spectra, biosig_interpolate
from transmission_spectra_simulator import TS_Simulator

from SEAS_Aux.calculation.interpolation import interpolate1d
import SEAS_Aux.calculation.astrophysics as calc 
import SEAS_Aux.cross_section.hapi as hp

import SEAS_Utils as utils
from SEAS_Utils.common_utils.constants import *
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.DIRs import *
from SEAS_Utils.common_utils.data_loader import two_column_file_loader,multi_column_file_loader, json_loader, molecule_cross_section_loader2
from SEAS_Utils.common_utils.data_saver import check_file_exist, check_path_exist
import SEAS_Utils.common_utils.db_management2 as dbm


class RS_Simulator(TS_Simulator):
    
    def __init__(self, user_input):
        TS_Simulator.__init__(self,user_input)
    
    
    def load_stellar_spectra(self):
        
        selection = self.user_input["Star"]["Spectra"]["selection"]
        filename  = self.user_input["Star"]["Spectra"]["filename"]
        
        orbit = utils.to_float(self.user_input["Planet"]["D_Planet"])
        R_P   = utils.to_float(self.user_input["Planet"]["R_Planet"])*R_Earth
        
        P_Surface = np.pi*R_P**2
        solar_constant_modifier = 1#(1./orbit)**2#*P_Surface*(0.01)
        
        data = multi_column_file_loader(os.path.join(Stellar_Spectra,selection,filename))
        x = np.array(data[0],dtype="float")
        y = np.array(data[1],dtype="float")*solar_constant_modifier

        normalized_flux = biosig_interpolate(10000./x,y,self.nu,"C")
        
        return normalized_flux
        
    def load_atmosphere_geometry_model(self, bio=False, CIA=False, Rayleigh=True):
        
        normalized_pressure         = self.normalized_pressure
        normalized_temperature      = self.normalized_temperature
        normalized_molecules        = self.normalized_molecules
        normalized_abundance        = self.normalized_abundance
        normalized_cross_section    = self.normalized_cross_section
        normalized_scale_height     = self.normalized_scale_height     

        if Rayleigh:
            normalized_rayleigh      = self.normalized_rayleigh     

    
        star_planet_observer_angle = 0
        incident_angle = star_planet_observer_angle*0.5
    
    
        # these parameters should be by wavelength/wavenumber ideally
        surface_scattering = 0.1
        surface_albedo = 0.9    # "glass like planet?"
        
        atmosphere_scattering = 0.1
        atmoshere_albedo = 0.01 # most like transmit through?
    
        

        Offset = float(self.user_input["Atmosphere_Effects"]["Base_Line"]["offset"])
        #Total_Transit_Signal = np.zeros(len(self.nu))+Offset
        
        
        # need to consider reflectance from each beam?
        # area?
        BeamTau = np.zeros(len(self.nu))
        for l,layer_height in enumerate(normalized_scale_height):
            
            pathl = layer_height/np.cos(incident_angle*np.pi/180.)
        
            ChunkTau = np.zeros(len(self.nu))
            for m, molecule in enumerate(normalized_molecules):
                
                molecular_ratio = normalized_abundance[l][m]
                number_density = (normalized_pressure[l]/(BoltK*normalized_temperature[l]))*molecular_ratio
                
                rayleigh = normalized_rayleigh[m]*molecular_ratio
                sigma = normalized_cross_section[m][l][l]
                
                
                effects = sigma+rayleigh
                
                # path length multiplied by 2 because the beam bounce back
                ChunkTau_Per_Molecule = number_density*(effects)*pathl*2*0.0001

                ChunkTau += ChunkTau_Per_Molecule     
    
            BeamTau += ChunkTau     
    
        Reflectance = calc.calc_reflectance(BeamTau,surface_albedo)
        
        # total receiving area?
        Raw_Total_Reflectance = Reflectance#*self.R_planet**2*(3./4)
    
    
        return Raw_Total_Reflectance
    
    
    
    
    
    
    