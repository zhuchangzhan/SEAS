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

Functions related to simulating observed spectra based on calculated theoratical spectra

for 0.8, need a full scale conversion of all list into dicts
instead of normalized_xxx, let's have a dict with pressure_layers as keys and relevent data as data

Takes in a simulated theoretical spectra and add observational effects


"""
import os
import sys
import numpy as np
import time
from scipy import interpolate, stats

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

#import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#ml = MultipleLocator(10)

import SEAS_Utils as utils
import SEAS_Aux.cross_section.hapi as hp
import SEAS_Main.observation_effects.noise as noise


class OS_Simulator():
    
    
    def __init__(self, user_input):
        
        self.user_input = user_input
     
    def add_noise(self, bin_mean_I, noise_type="gaussian"):
        
        error_scale = utils.to_float(self.user_input["Observation_Effects"]["Noise"]["error_scale"])
        
        data_length = len(bin_mean_I)
        
        if noise_type == "gaussian":
            Noise = noise.Gaussian_Noise(error_scale, data_length)
            noise_added = Noise.get_noise()
            
        elif noise_type == "poisson":
            Noise = noise.Poisson_Noise(error_scale, data_length) 
            noise_added = Noise.get_noise()

        bin_mean_error_I = bin_mean_I*(1.0+noise_added)
        bin_mean_error_bar = error_scale*bin_mean_error_I[1]  
        
        
        return bin_mean_error_I, bin_mean_error_bar
        

    def add_background_stars(self):
        pass

    def calculate_convolve(self, nu, trans, record=False):
        #if self.user_input["Observation"]["Convolve"] == "true":
        amount = utils.to_float(self.user_input["Observation_Effects"]["Convolve"]["convolve_amount"])
        nu,Transit_Signal,i1,i2,slit = hp.convolveSpectrum(nu,trans,SlitFunction=hp.SLIT_RECTANGULAR,Resolution=amount,AF_wing=20.0)
        
        
        if record:
            return nu,Transit_Signal,max(Transit_Signal)
        
        return nu,Transit_Signal

    def calculate_bin(self, x, signal):
      
        
        Bins        = utils.to_float(self.user_input["Observation_Effects"]["Bin"]["bin_number"])
        Method      = self.user_input["Observation_Effects"]["Bin"]["method"]
        
        bin_mean_I, bin_edges, binnumber = stats.binned_statistic(x, signal, statistic=Method, bins=Bins)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2
        
        
        return bin_centers, bin_mean_I
    
    def spectra_to_magnitude(self):
        """
        convert simulated spectra to magnitude
        This will need a simulated stellar spectra (a black body curve)?
        """
        
        pass
    
    def number_of_photon(self):
        """
        number of photons expected per bin
        """
        pass
    
    
    def telescope_response_function(self):
        pass
    
    
    def telescope_jitter(self):
        pass




