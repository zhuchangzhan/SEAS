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
Cross Section Calculator

Molecular Line List from HITRAN are read in and cross section are calculated using the hapi.py from HITRAN

"""

import SEAS_Aux.cross_section.hapi as hp
from SEAS_Utils import to_float,to_int
import lines2xsec as l2x


class cross_section_calculator():


    def __init__(self,d_path,r_path,molecule,component,numin,numax,step=0.1,P=1,T=300,
                 gamma="gamma_self",cross=True, direct=False, multi=False):
        """
        
        You idiot why no P unit?? This is probably 1 Pa
        
        n: molecule name
        m: molecular isotope
        i: molecule abundance
        
        r_path is not used. This is reserved for saving the individual cross sections
        no need to implement at this point
        
        """

        self.d_path = d_path
        self.r_path = r_path
        
        self.molecule = molecule
        
        self.n = component[0]
        self.m = component[1]
        self.i = component[2]

        self.P = P
        self.T = T

        self.numin = numin
        self.numax = numax   
        self.step = step
        
        self.gamma = gamma
        self.cross = cross
        
        self.direct = direct
        self.multi = multi
        
        #self.fetch_data()
    
    
    def fetch_data(self):
        """
        data input for hapi calculator
        """
        
        hp.db_begin(self.d_path)
        hp.select(self.n,ParameterNames=("nu","sw"), Conditions=("between","nu",to_float(self.numin),to_float(self.numax)))
    
    
    
    def hapi_calculator(self):
        """
        This calculation method still have speed issues
        Also issues with the select function in hapi not working as intended
        
        Not used in 0.7, only as placeholder
        """
        
        nu, coef = hp.absorptionCoefficient_Voigt(((self.n,self.m,self.i),),
                                                  self.molecule, 
                                                  OmegaStep=self.step,
                                                  HITRAN_units=self.cross,
                                                  GammaL=self.gamma, 
                                                  Environment={'p':to_float(self.P),'T':to_float(self.T)})
    
        return nu, coef
    

    def read_data(self):
        """
        data input for personal_calculator
        """
    
        molecule_data = l2x.read_data(self.d_path,self.molecule,self.numin,self.numax,direct=self.direct,multi=self.multi)
        
        return  molecule_data
    
    
    def personal_calculator(self, data = None):
        """
        This is an inherited calculation method from the hapi calculator and offers a bandage patch
        approach solution to problems mentioned above. This is the line_list_to_cross_section_calculation
        code used in version 0.4-0.6 and now re-wrapped under 0.7 guide lines.
        
        There are also problems with this calculation method. Will need a complete look in 0.8
            can we remove the component? 
            if possible, implement ctop 
            a working, speedy hapi calculator?
            why do we need self.numin and numax again for calculation when we've already constrained it in read_data?...
            or... should we read in a bit more by default?
            
        How is this compatible with other line list databases?
        
        """
        if data == None:
            molecule_data = self.read_data()
        else:
            molecule_data = data
            
        component = [self.n, self.m, self.i]
        
        nu, coef = l2x.absorption_Voigt_calculation(molecule_data, component, 
                                                    self.gamma, self.P, self.T, 
                                                    self.numin, self.numax, self.step)

        return nu, coef












