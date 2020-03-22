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
This contains the astrophysics equations of the atmosphere simulation

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.Common_Utils.constants import *
from SEAS_Main.Physics.molecular_weight import calculate_mw

def blackbody_lam(wav, T):
    """ Blackbody as a function of wavelength (m) and temperature (K).
    returns units of erg/s/cm^2/cm/Steradian
    """
    a = 2*HPlanck*CLight**2
    b = HPlanck*CLight/(wav*BoltK*T)
    intensity = a/((wav**5)*(np.exp(b)-1.0))
    return intensity

def blackbody_nu(wn,T): # nu is frequency (Hz) , T is temperature (K)
    """ 
    Blackbody as a function of wavenumber (cm^-1) and temperature (K).
    """
    nu = wn*3e10 # convert wavenumber to frequency
    return (2*h*nu**3/c**2)*(1/(np.exp((h*nu)/(k*T))-1))

def planck(wav,T):
    """
    calculate the planck's blackbody intensity
    redundant with blackbody_lam?
    """
    return (2*HPlanck*c/wav**3)*(1/(np.exp((HPlanck*c)/(wav*BoltK*T))-1))

def calc_transmittance(tau):
    """
    calculate transmission from the equation e^(-tau)
    """
    return np.exp(-tau)

def calc_reflectance(tau, albedo):
    
    return np.exp(-tau)*albedo

def calc_H(Temperature, MeanMolWeight, SurfaceG, r=2):
    """
    calculate atmospheric scale height
    """
    
    return round((BoltK*Temperature)/(MeanMolWeight*SurfaceG),r)

def get_MolWeight(molecules):
    """
    get the molecular weight for molecules
    """
    mw = np.zeros(len(molecules))
    for i,molecule in enumerate(molecules):
        mw[i] = calculate_mw(molecule)
    return mw
        
    
    

def calc_MeanMolWeight(molecules, molecular_weight):
    """
    calculate the mean molecular weight
    """
    
    total = 0
    for i in range(len(molecules)):
        total+= molecules[i]*molecular_weight[i]
    return total
    
def calc_SurfaceG(M_planet,R_planet):
    """
    calculate the surface gravity
    """
    return G*M_planet/R_planet**2

def calc_EqTemperature(distance):
    """
    calculate the equilibrium temperature at the top of the atmosphere.
    
    Not implemented yet. Will look back when considering non-isothermal atmospheres
    """
    return 0



def calc_rayleigh(molecule, nu_range):
    """
    calculate rayleigh scattering
    is N a constant? really? 
    
    """
    #N = N/10**6
    
    N = 2.547*10**19
    

    front = 24*np.pi**3*nu_range**4/N**2
    
    if molecule == "CH4":
        # 46662+-13, 4.02+-0.03. how to implement this?
        nvmo = 46662 + 4.02*10**-6*nu_range**2
        nv = nvmo/10**8 +1
        Fk = 1.0
    
    elif molecule == "CO":
        #22851+-200, 0.456+-0.01, 71427+-200
        nvmo = 22851+(0.456*10**12)/(71427**2-nu_range**2)
        nv = nvmo/10**8 +1
        Fk = 1.016    
    
    elif molecule == "CO2":
        return np.zeros(len(nu_range))
        a = 5799.25/(128908.9**2-nu_range**2)
        b = 120.05/(89223.8**2-nu_range**2)
        c = 5.3334/(5037.5**2-nu_range**2)
        d = 4.3244/(67837.7**2-nu_range**2)
        e = 0.1218145*10**-4/(2418.136**2-nu_range**2)
        nvmo = 1.1427*10**6*(a+b+c+d+e)
        nv = nvmo +1
        # 1.1364+-0.0005, 25.3+-1.5
        Fk = 1.1364+25.3*10**-12*nu_range**2   
        
    elif molecule == "N2O":
        # 46890+-85, 4.12+-0.2
        nvmo = 46890+4.12*10**-6*nu_range**2
        nv = nvmo/10**8 +1
        Fk = 1.225    

    elif molecule == "N2":
        # what to do if nu<4860 and nu>39370?
        
        numin = nu_range[0]
        numax = nu_range[-1]
        
        
        if numax <= 21360:
            nvmo = 6498.2 + (307.43305*10**12)/(14.4*10**9-nu_range**2) 
        elif numin >= 21360:
            nvmo = 5677.465 + (318.81874*10**12)/(14.4*10**9-nu_range**2)
        else:
            nu_range_low = []
            nu_range_high = []
            for i in nu_range:
                if i < 21360:
                    nu_range_low.append(i)
                else:
                    nu_range_high.append(i)
            
            nu_range_low = np.array(nu_range_low)
            nu_range_high = np.array(nu_range_high)
            
            nvmol = 6498.2 + (307.43305*10**12)/(14.4*10**9-nu_range_low**2)
            nvmoh = 5677.465 + (318.81874*10**12)/(14.4*10**9-nu_range_high**2)
            nvmo = np.concatenate([nvmol,nvmoh])    
            
        
        Fk = 1.034+3.17*10**-12*nu_range**2    
        nv = nvmo/10**8 +1
           
        
    elif molecule == "O2":
        nvmo = 20564.8+(24.80899)/(4.09*10**9-nu_range**2)
        nv = nvmo/10**8 +1
        Fk = 1.096+1.385*10**-11*nu_range**2+1.448*10**-20*nu_range**2    
    
    elif molecule == "H2":
        # not implemented yet
        lambd = 10000./nu_range*10**-4
        sigma = 8.49*10**-45/lambd**4
        
        return sigma
        
    else:
        return np.zeros(len(nu_range))

    
    
    
    middle = ((nv**2-1)/(nv**2+2))**2    
    rayleigh = front*middle*Fk
    
    return rayleigh 


def calc_CIA(molecule, nu):
    """
    Calculate Collision Induced Opacities
    """
    
    return np.zeros(len(nu))

def calc_cloud(height, nu):
    """
    Calculate Cloud Model
    """
    
    return np.zeros(len(nu))
    
def calc_cloud_number_density(air_number_density = 1.225e-3, # g/cm^3
                              particle_mixing_ratio = 4.62e-6, #in abs abundance, unitless
                              particle_density = 4.09, # g/cm^3
                              particle_radius = 1e-4): # cm
    
    # assuming spherical particles
    unit_particle_mass = 4/3*np.pi*particle_density*particle_radius**3 #g/particle
    
    particle_vapor_density = air_number_density*particle_mixing_ratio #g/cm^3

    particles_number_density = particle_vapor_density/unit_particle_mass #particle/cm^3
    
    return particles_number_density # number of particles/cm^3

def calculate_air_density(P,T,mean_air=28):
    
    # stuff these two into constants
    k = 1.38e-23
    R = 6.022e23
    return P*mean_air/(k*T*R)
        
    
    