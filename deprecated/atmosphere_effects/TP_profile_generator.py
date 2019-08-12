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
This is a logistic module that generates the Temperature Pressure Profiles
This module does not handle the physics and chemistry that goes into 
determining the actual Temperature Pressure Profiles of the atmosphere

The exact purpose of this generator remains unclear until working with
atmosphere physics and chemistry simulations...
For now it's a placeholder for future code

"""
import os
import sys
import numpy as np

from SEAS_Utils import to_float
from SEAS_Utils.common_utils.data_saver import check_file_exist, check_path_exist

class temperature_pressure_profile_generator():
    
    
    def __init__(self,
                 data,
                 pressures = [100000,10000,1000,100,10,1,0.1,0.01,0.001,0.0001,0.00001],
                 path = "../../input/atmosphere_data/TP_Profile",
                 name = "Temp.txt",
                 overwrite = False
                 ):
        
        
        self.data      = data
        self.pressures = pressures
        self.path      = path
        self.name      = name
        self.overwrite = overwrite
    
    def load(self):
        """
        reserve for biochemistry study use
        """
        pass
    
    
    
    def generate(self):
        
        T_Surface = to_float(self.data["Surface_Temperature"])
        Type      = self.data["Type"]
        
        T = []
        if Type == "isothermal":
            for _ in range(len(self.pressures)):
                T.append(T_Surface)
        elif Type == "variable":
            for _ in range(len(self.pressures)):
                T.append(T_Surface)
            
            for change in self.data.keys()[2:]:
                T_Start   = to_float(self.data[change]["Start_Temperature"])
                T_End     = to_float(self.data[change]["End_Temperature"])
                P_Start   = to_float(self.data[change]["Start_Pressure"])
                P_End     = to_float(self.data[change]["End_Pressure"])
                
                for i,pres in enumerate(self.pressures):
                    if pres >= P_Start:
                        pass
                    elif pres <= P_End:
                        T[i] = T_End
                    else:
                        T[i] = T_Start + (T_End-T_Start)*np.log10(P_Start/pres)/np.log10(P_Start/P_End)
                        
                    
        self.T = T
        
        
        print "TP_Profile generated!"
        return T
        



    def save(self):
        
        check_path_exist(self.path)
        save_path = os.path.join(self.path,self.name)
        check_file_exist(save_path)
        
        print self.T
        print self.pressures
        
        with open(save_path,"w") as file:
            for i,t in enumerate(self.T):
                file.write(" ".join([str(self.pressures[i]),str(t)]))
                if i == len(self.T)-1:
                    break
                file.write("\n")

        print "TP_Profile saved to %s"%save_path


