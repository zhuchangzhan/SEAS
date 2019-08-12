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

GUI for using the SEAS

for toggable, dynamic plotting of SEAS simulations


simulation presets & description

develop more on the plotter
have multiple plot displayed.

Not a Spectral analyzer but a "Result Analyzer"





"""


import os
import sys
import shutil
import numpy as np
from Tkinter import *

import PIL
from PIL import Image, ImageTk

import matplotlib
from matplotlib.figure import Figure
matplotlib.use('TkAgg')
from matplotlib import pylab as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze

import SEAS_Utils as utils
import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.timer import simple_timer




class GUI(Frame):
    
    def __init__(self, master):
        Frame.__init__(self,master)



class SEAS_Main_GUI(GUI):
    
    def __init__(self, master):
        GUI.__init__(self, master) 

        


    def main_frame(self):
        """



        molecule selection

        


        """
        
        pass 






