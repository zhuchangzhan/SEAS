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
Module for loading and utilizing CIA data

This reader currently works for a subset of the simulation
currently omitting the eq and norm versions

"""
import os
import sys
import itertools
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.DIRs import HITRAN_CIA



def find_data():
    data_file = []
    for file in os.listdir(HITRAN_CIA):
        if "cia" in file and len(file.replace("_","-").split("-"))==3:
            data_file.append(file)
    return data_file


def select_molecule_cia(molecule_list,user_defined="-1#@%"):
    """
    various ways of determining which cia files to load
    """

    combination = ["-".join([x,x]) for x in molecule_list] + \
                  ["-".join(x) for x in itertools.combinations(molecule_list,2)] + \
                  ["-".join([x[1],x[0]]) for x in itertools.combinations(molecule_list,2)] + \
                  [user_defined]

    datafile = find_data()
    
    result = []
    for i in combination:
        for j in datafile:
            if i in j:
                result.append(j)
    return result



class HITRAN_CIA_data_processor():
    
    def __init__(self,filepath,filename):
        
        self.filepath = filepath
        self.filename = filename
        
    def load(self): 
        
        with open(os.path.join(self.filepath,self.filename)) as f:
    
            self.temperature, self.nu, self.xsec = [],[],[]

            try:# files with single datapoint length
                data = f.read().split("\n")
                if data[-1] == "":
                    data = data[:-1]
                
                data_point = int(data[0][40:47].strip())
                data_grid = np.reshape(np.array(data),(-1,data_point+1))
                for i in data_grid:
                    self.temperature.append(float(i[0][47:54].strip()))
                    n,x = [list(x) for x  in zip(*[x.split() for x in i[1:]])]
                    self.xsec.append(np.array(x,dtype="float"))

                self.nu = n 
            
                """
                formula      = header[:20].strip()
                numin        = float(header[20:30].strip())
                numax        = float(header[30:40].strip())
                data_point   = int(header[40:47].strip())
                temperature  = float(header[47:54].strip())
                maximum_xsec = float(header[54:64].strip())
                resolution   = float(header[64:70].strip())
                comments     = header[76:97].strip()
                reference    = header[97:].strip()
                """
            
            except: # files with different datapoint length. Not sure what to do now
                pass            
            
            return self.temperature, self.nu, self.xsec

            
    
    
    