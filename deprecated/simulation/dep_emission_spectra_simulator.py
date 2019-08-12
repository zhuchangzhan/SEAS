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

Simple Emission Spectra Simulator 

"""
import os
import sys
import numpy as np
import time
from scipy import interpolate
#import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#ml = MultipleLocator(10)

import SEAS_Main.atmosphere_geometry
import SEAS_Main.atmosphere_property
from SEAS_Main.atmosphere_effects.biosig_molecule import *
from SEAS_Main.atmosphere_effects.cloud import *

from SEAS_Aux.calculation.interpolation import interpolate1d
import SEAS_Aux.calculation.astrophysics as calc 

# DIRs should actually be a toggable param in the cfg file
import SEAS_Utils.common_utils as utils
from SEAS_Utils.common_utils.DIRs import *
from SEAS_Utils.common_utils.constants import *
from SEAS_Utils.common_utils.data_loader import *
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.data_saver import check_file_exist, check_path_exist
import SEAS_Utils.common_utils.db_management2 as dbm

from scipy.constants import h,k,c

class ES_Simulator():
    
    def __init__(self, user_input):
        self.user_input = user_input

    def blackbody_lam(self, lam, T):
        """ Blackbody as a function of wavelength (um) and temperature (K).
    
        returns units of erg/s/cm^2/cm/Steradian
        """
    
        lam = 1e-6 * lam # convert to metres
        return 2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))









