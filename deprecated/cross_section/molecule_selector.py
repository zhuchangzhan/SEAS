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
Due to our data coming from multiple sources, we need a way to automatically   
select which molecule gets used in the simulation.

"""

from SEAS_Utils.common_utils.DIRs import *
from SEAS_Utils.common_utils.data_saver import check_path_exist, check_file_exist




def main_molecule_selector(molecule, 
                           preference_order=["Exomol","HITRAN_Lines","HITRAN_Xsec","NIST"],
                           auxillary=True):

    for preference in preference_order:
        
        if preference == "Exomol":
            
            
            
            
        elif preference == "HITRAN_Lines":
            
        elif preference == "HITRAN_Xsec":
            
        elif preference == "NIST":
            
        
        
        
        
        
        print preference



    #return selection













