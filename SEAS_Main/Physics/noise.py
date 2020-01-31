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
Noise Module. Handles all noise associated works

Field convention is to convert the spectra into some sort of flux and then photon count
Then use sqrt of photon count as noise base?

or a defined noise floor?

or ... stuff. 

What exactly is JWST Pandexo simulation doing anyway?

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Utils.External_Utils.hapi as hp 

def convolve_spectra():
    pass

def simple_noise():
    
    
    pass

class Noise():
    
    def __init__(self):
        pass
    
    def get_noise(self):
        pass
    
    


class Photon_Noise():
    
    def __init__(self, user_input):
        
        self.user_input = user_input

    def blackbody_lam(self, wav, T):
        """ Blackbody as a function of wavelength (m) and temperature (K).
        returns units of erg/s/cm^2/cm/Steradian
        """
        a = 2*HPlanck*CLight**2
        b = HPlanck*CLight/(wav*BoltK*T)
        intensity = a/((wav**5)*(np.exp(b)-1.0))
        return intensity

    def calculate_noise(self, D_atmosphere):
    
        R_Star       = self.user_input["Star"]["R_Star"]
        R_planet     = self.user_input["Planet"]["R_Planet"] 
        
        if R_Star < 1000:
            R_Star *= R_Sun
        
        if R_planet < 1000:
            R_planet *= R_Earth
        
        
        #D_atmosphere = np.ones(len(D_atmosphere))*max(D_atmosphere)
        
        R_obs        = self.user_input["Telescope"]["Aperture"]
        #R_atmosphere = self.user_input["Planet"]["R_Atmosphere"]
        Distance     = self.user_input["Telescope"]["Distance"]
        Duration     = self.user_input["Telescope"]["Duration"]
        Quantum      = self.user_input["Telescope"]["Quantum_Efficiency"]
        T_Star       = self.user_input["Star"]["T"]
        Noise_M      = self.user_input["Observation_Effects"]["Noise"]["Multiplier"]
        
        
        
        # calculate number of photons 
        B_Body = self.blackbody_lam(self.bin_centers*10**-6, T_Star)
        Bin_width = self.bin_width*10**-6
        A_Star = np.pi*R_Star**2
        Psi_Tele = np.pi*R_obs**2/Distance**2
        E_Total = B_Body*Bin_width*A_Star*Psi_Tele*Duration
        num_photon = (E_Total*self.bin_centers*10**-6)/(HPlanck*c)*Quantum
        
        #signal = (2*R_planet*D_atmosphere)/R_Star**2*num_photon
        
        
        #R_planet = R_planet+D_atmosphere
        
        signal = R_planet**2/R_Star**2*num_photon
        photon_noise = Noise_M*np.sqrt(num_photon)
        SNR = signal/photon_noise

                
        return signal, photon_noise, SNR
    
    def determine_bin(self):
    
        lambda_init    = float(self.user_input["Telescope"]["min_wavelength"])
        lambda_max     = float(self.user_input["Telescope"]["max_wavelength"])
        bin_width_init = float(self.user_input["Telescope"]["Binning"]["bin_width"])
        bin_exponent   = float(self.user_input["Telescope"]["Binning"]["bin_exponent"])
        
        bin_edges,bin_width,bin_centers = [],[],[]
        
        
        
        i=0
        lambda_current = lambda_init
        bin_edges.append(lambda_current)
        while True:
            
            new_bin_width = bin_width_init*(lambda_current/lambda_init)**(bin_exponent)
            
            lambda_center = lambda_current+0.5*new_bin_width
            lambda_current += new_bin_width
            
            if lambda_center > lambda_max:
                break
            bin_edges.append(lambda_current)
            bin_centers.append(lambda_center)
            bin_width.append(new_bin_width)
            
            
            i+=1
        
        self.bin_edges = np.array(bin_edges)
        self.bin_width = np.array(bin_width)
        self.bin_centers = np.array(bin_centers)
        
        
        return self.bin_edges, self.bin_width, self.bin_centers

    def calculate_convolve(self, nu, trans):
        #if self.user_input["Observation"]["Convolve"] == "true":
        amount = int(self.user_input["Spectra"]["convolve_amount"])
        nu,Transit_Signal,i1,i2,slit = hp.convolveSpectrum(nu,trans,SlitFunction=hp.SLIT_RECTANGULAR,Resolution=10,AF_wing=20.0)
        
        return nu,Transit_Signal

    def calculate_bin(self, x, signal):
      
        
        Bins        = utils.to_float(self.user_input["Observation_Effects"]["Bin"]["bin_number"])
        Method      = self.user_input["Observation_Effects"]["Bin"]["method"]
        
        bin_mean_I, bin_edges, binnumber = stats.binned_statistic(x, signal, statistic=Method, bins=Bins)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2
        
        
        return bin_centers, bin_mean_I




class Poisson_Noise(Noise):
    
    def __init__(self, poisson, shape):
        
        poissonNoise = np.random.poisson(poisson, shape).astype(float)
    
        return poissonNoise
    
class Shot_Noise(Poisson_Noise):
    """
    shot noise is most frequently observed with small currents or low light intensities that have been amplified.
    
    SNR = sqrt(N), where N is average number of event
    
    """    
    def __init__(self):
        pass
        
class Gaussian_Noise(Noise):

    def __init__(self, multiplier=1, length=10):
    
        self.length = length
        self.multiplier = multiplier
    
    def get_noise(self):
        return np.random.randn(self.length)*self.multiplier
        
        """
        mu, sigma = 8560, 20 # mean and standard deviation
        s = np.random.normal(mu, sigma, 1000)
        
        import matplotlib.pyplot as plt
        count, bins, ignored = plt.hist(s, 250, normed=True)
        print count, bins
        
        
        func = 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) )
        
        plt.plot(bins, func, linewidth=2, color='r')
        plt.show()
        """
        
class Uniform_Noise(Noise):

    def __init__(self):
        
        
        return np.random.random()
    
class Laplace_Noise(Noise):

    def __init__(self):
        pass
    
class Lorentz_Noise(Noise):

    def __init__(self):
        pass

class Perlin_Noise(Noise):
    """
    reference https://pypi.python.org/pypi/noise/
    """
    
    def __init__(self):
        pass

class Telescope_Noise(Noise):

    def __init__(self):
        pass


    def add_jitter(self):
        pass





    
    