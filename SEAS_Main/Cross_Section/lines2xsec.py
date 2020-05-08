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
A Re-worked version of the old line_list_to_cross_section_calculation used in 0.4-0.6 
to generated molecular cross sections. 

This is still a bandage patch and does not solve some of the core problems with the fix

Will need a complete new look to check for accuracy in version 0.8

"""

import os
import hapi as hp
import numpy as np
from bisect import bisect
from numpy import complex128,int64,float64,exp

SI = False

# default values for intensity threshold
DefaultIntensityThreshold = 0. # cm*molec
# default value for omega wing
DefaultOmegaWing = None
# default value for omega wing in halfwidths (from center)
DefaultOmegaWingHW = 50. # cm-1    HOTW default

# Defining precision
__ComplexType__ = complex128
__IntegerType__ = int64
__FloatType__ = float64

if SI:
    cBolts = 1.3806503e-23
    cc = 2.99792458e8
else: #CGS
    cBolts = 1.380648813E-16 # erg/K, CGS
    cc = 2.99792458e10 # cm/s, CGS

# No Bias
cMassMol = 1.66053873e-27
cZero = 0.

# Hapi dependencies
# this is too much work to isolate and handle
PROFILE_VOIGT = hp.PROFILE_VOIGT   
PYTIPS=hp.PYTIPS


# Volume concentration of all gas molecules at the pressure p and temperature T
def volumeConcentration(p,T):
    if SI:
        return (p/9.869233e-6)/(cBolts*T) # SI
    else:
        return (p/9.869233e-7)/(cBolts*T) # CGS

# temperature dependence for intencities (HITRAN)
def EnvironmentDependency_Intensity(LineIntensityRef,T,Tref,SigmaT,SigmaTref,
                                    LowerStateEnergy,LineCenter):
    const = __FloatType__(1.4388028496642257)  # check this unit
    ch = exp(-const*LowerStateEnergy/T)*(1-exp(-const*LineCenter/T))
    zn = exp(-const*LowerStateEnergy/Tref)*(1-exp(-const*LineCenter/Tref))
    LineIntensity = LineIntensityRef*SigmaTref/SigmaT*ch/zn
    return LineIntensity

def EnvironmentDependency_Gamma0(Gamma0_ref,T,Tref,p,pref,TempRatioPower):
    return Gamma0_ref*p/pref*(Tref/T)**TempRatioPower


def read_data(path,molecule,numin,numax,imin=-1,ctop=-1,direct=False,multi=False):
    """
    numin: minimum wavenumber
    numax: maximum wavenumber
    imax:  minimum intensity
    ctop:  Top intensity  
    
    Reading data from locally stored line list files 
    Can be constrained in wavenumber
    Returns a dictionary containing data from the line list file
    """

    if multi:
        
        M, I, wavenumber, intensity,gamma_air,gamma_self, elower, n_air,delta_air = [],[],[],[],[],[],[],[],[]
            
        for filename in path:
            
            front,back = filename.split("-")
            
            filemin = float(front.split("_")[-1])
            filemax = float(back.split("_")[0])
            
            if filemin > numax or filemax < numin:
                continue
            
            
            #print "Reading Data From %s"%(filename)
            
            
            f = open(filename).read().split("\n")
            
            
            for i in range(len(f)-1):
                
                wn = float(f[i][3:15])
                
                if wn < numin or wn >= numax:
                    continue
            
                M.append(         int(f[i][0:2]))
                I.append(         int(f[i][2]))
                wavenumber.append(float(f[i][3:15]))
                intensity.append( float(f[i][16:26]))
                gamma_air.append( float(f[i][35:40]))
                gamma_self.append(float(f[i][41:45]))
                elower.append(    float(f[i][46:56]))
                n_air.append(     float(f[i][56:59]))
                delta_air.append( float(f[i][60:68]))
                
        data = {"M":M, "I":I, "wavenumber":wavenumber, "intensity":intensity, 
                  "gamma_air":gamma_air, "gamma_self":gamma_self, "elower":elower,
                  "n_air":n_air, "delta_air":delta_air}
        
        return data
    
    else:
        
        if direct:
            print("Reading Data From %s"%(path))
            f = open(path).read().split("\n")
            
        else:
            #print "Reading Data From %s/%s.data"%(path,molecule)
            f = open("%s/%s.data"%(path,molecule)).read().split("\n")
        
        M, I, wavenumber, intensity,gamma_air,gamma_self, elower, n_air,delta_air = [],[],[],[],[],[],[],[],[]
        
        appending = False
        for i in range(len(f)-1):
            
            wn = float(f[i][3:15])
            
            if wn >= numin:
                appending = True
            if wn >= numax:
                break
            
            M.append(         int(f[i][0:2]))
            I.append(         int(f[i][2]))
            wavenumber.append(float(f[i][3:15]))
            intensity.append( float(f[i][16:26]))
            gamma_air.append( float(f[i][35:40]))
            gamma_self.append(float(f[i][41:45]))
            elower.append(    float(f[i][46:56]))
            n_air.append(     float(f[i][56:59]))
            delta_air.append( float(f[i][60:68]))
            
        data = {"M":M, "I":I, "wavenumber":wavenumber, "intensity":intensity, 
                  "gamma_air":gamma_air, "gamma_self":gamma_self, "elower":elower,
                  "n_air":n_air, "delta_air":delta_air}
        
        return data


def absorption_Voigt_calculation(molecule_data, component, gamma_name, P, T, numin, numax, step=0.1, cross_section=True, OmegaWing=DefaultOmegaWing,
                 IntensityThreshold=DefaultIntensityThreshold, OmegaWingHW=DefaultOmegaWingHW, partitionFunction=PYTIPS):
    """
    Considering both Doppler Broadening and Pressure Broadening
    
    calculating cross section for each molecule
    
    P: Pressure in Pa, although need to look into P. 1atm should be 101300. For simplicity we're using 100000... for now.
    
    """
    
    # Reference Conditions
    Tref = 296.
    Pref = 1.
    P = P/100000.0
    
    
    # Determining the length of the wavenumber array
    nu = molecule_data["wavenumber"]
    
    
    Omega_min = float(numin)
    Omega_max = float(numax)
    OmegaStep = step   
    OmegaCount = (Omega_max-Omega_min)/OmegaStep+1
    Omegas = np.linspace(Omega_min,Omega_max,OmegaCount)[:-1]
    Xsect  = np.zeros(len(Omegas)) 
    
    
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    

    M = molecule_data["M"][0]
    I = molecule_data["I"][0]
    if len(component) >=3:
        N = component[2]
    else:
        N = hp.ISO[(M,I)][hp.ISO_INDEX["abundance"]]
    
    ABUNDANCES[(M,I)] = N
    NATURAL_ABUNDANCES[(M,I)] = hp.ISO[(M,I)][hp.ISO_INDEX["abundance"]]    

    
    factor = hp.volumeConcentration(P,T)

    EnvDependences = lambda ENV, LINE:{}
    Env = {}
    Env["T"] = float(T)
    Env["P"] = float(P)
    Env["Tref"] = Tref
    Env["Pref"] = Pref    
    

        
    nline = len(molecule_data["wavenumber"])
    MolNumberDB   = int(molecule_data["M"][0])
    IsoNumberDB   = int(molecule_data["I"][0])
    
    if SI:
        m = hp.ISO[(MolNumberDB,IsoNumberDB)][hp.ISO_INDEX['mass']]*cMassMol* cMassMol 
    else:
        m = hp.ISO[(MolNumberDB,IsoNumberDB)][hp.ISO_INDEX['mass']]*cMassMol* 1000.
    
    parnames = molecule_data.keys()
    
    # alot of calculations can be omitted here given we don't need to calculate the entire thing
    for RowID in range(nline):
        
        Line = {parname:molecule_data[parname][RowID] for parname in parnames}
        CustomEnvDependences = EnvDependences(Env,Line)
        
        if (MolNumberDB,IsoNumberDB) not in ABUNDANCES: continue
        
        SigmaT = partitionFunction(MolNumberDB,IsoNumberDB,T)
        SigmaTref = partitionFunction(MolNumberDB,IsoNumberDB,Tref)   
         
        LineCenterDB = molecule_data['wavenumber'][RowID]
        
        
        LineIntensityDB = molecule_data['intensity'][RowID]
        LowerStateEnergyDB = molecule_data['elower'][RowID]

        if 'intensity' in CustomEnvDependences:
            LineIntensity = CustomEnvDependences['intensity']
        else:
            LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB,T,Tref,SigmaT,SigmaTref,
                                                            LowerStateEnergyDB,LineCenterDB)
        
        #   FILTER by LineIntensity: compare it with IntencityThreshold
        if LineIntensity < IntensityThreshold: continue

        # Doppler Broadening
        GammaD = np.sqrt(2*cBolts*T*np.log(2)/m/cc**2)*LineCenterDB
        
        # Pressure Broadening
        Gamma0 = 0.
        Shift0 = 0.

        # default is air, not considering multiple diluent atm
        abun = 1
        try:
            Gamma0DB =  molecule_data[gamma_name][RowID]
        except:
            Gamma0DB = 0.0
        
        
        if gamma_name == "gamma_air":
            n_name = "n_air"
            d_name = "delta_air"
            dp_name= "deltap_air"
        elif gamma_name == "gamma_self":
            n_name = "n_self"
            d_name = "delta_self"
            dp_name= "deltap_self"
        else:
            print("unknown gamma")
            return
        
        try:
            TempRatioPowerDB = molecule_data[n_name][RowID]
            if n_name == "n_self" and TempRatioPowerDB == 0.:
                TempRatioPowerDB = molecule_data["n_air"][RowID]
        except:
            TempRatioPowerDB = molecule_data["n_air"][RowID]
        
        Gamma0 += abun*CustomEnvDependences.get(gamma_name,
                  EnvironmentDependency_Gamma0(Gamma0DB,T,Tref,P,Pref,TempRatioPowerDB))

        try:
            Shift0DB = molecule_data[d_name][RowID]
        except:
            Shift0DB = 0.0
        
        try:
            deltap = molecule_data[dp_name][RowID]
        except:
            deltap = 0.0
        
        # pressure broadening for self is 0?
        Shift0 += abun*CustomEnvDependences.get(d_name, # default ->
                  ((Shift0DB + deltap*(T-Tref))*P/Pref))
        #Shift0 = 0
    
        OmegaWingF = max(OmegaWing,OmegaWingHW*Gamma0,OmegaWingHW*GammaD)
        BoundIndexLower = bisect(Omegas,LineCenterDB-OmegaWingF)
        BoundIndexUpper = bisect(Omegas,LineCenterDB+OmegaWingF)
        
        lineshape_vals = PROFILE_VOIGT(LineCenterDB+Shift0,GammaD,Gamma0,Omegas[BoundIndexLower:BoundIndexUpper])[0]
        
        
        # absorption coefficient calculation
        
        
        if cross_section:
            Xsect[BoundIndexLower:BoundIndexUpper] += ABUNDANCES[(MolNumberDB,IsoNumberDB)]/NATURAL_ABUNDANCES[(MolNumberDB,IsoNumberDB)] * \
                                                      LineIntensity * lineshape_vals
        else:
            Xsect[BoundIndexLower:BoundIndexUpper] += factor / NATURAL_ABUNDANCES[(MolNumberDB,IsoNumberDB)] * \
                                                      ABUNDANCES[(MolNumberDB,IsoNumberDB)] * \
                                                      LineIntensity * lineshape_vals

        
    return Omegas, Xsect






