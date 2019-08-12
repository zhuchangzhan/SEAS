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


atmosphere window/filter nomenclature?

"""
import os
import sys
import numpy as np
import time
from scipy import interpolate


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

#import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#ml = MultipleLocator(10)

import SEAS_Aux.cross_section.hapi as hp
import SEAS_Main.observation_effects.noise as noise



class Spectra_Analyzer():
    
    
    def __init__(self, user_input):
        
        self.user_input = user_input
    
    def new_spectra_window(self,nu,coef,ratio=0.3,span=100,result="nu",thres=False):
        """
        window should be evaluated at different regions
        
        < 1 um cut out cus jwst won't need it
        1-5um [10000-2000]
        5-15um [2000-666]
        15-25um [666-400]
        """
        windows = []
        vis,near,mid,far = [],[],[],[]
        VIS,NEAR,MID,FAR = [],[],[],[]
        
        for i,n in enumerate(nu):
            if n > 10000:
                vis.append(n)
                VIS.append(coef[i])
            elif n > 2000 and n <= 10000:
                near.append(n)
                NEAR.append(coef[i])
            elif n > 666 and n <= 2000:
                mid.append(n)
                MID.append(coef[i])
            elif n > 400 and n <= 666:
                far.append(n)
                FAR.append(coef[i])
        
        if thres == False:
            near_window = self.find_window(near,NEAR,ratio,span,result,thres)
            mid_window  = self.find_window(mid,MID,ratio,span,result,thres)
            far_window  = self.find_window(far,FAR,ratio,span,result,thres)
            
            for win in near_window:
                windows.append(win)
            for win in mid_window:
                windows.append(win)
            for win in far_window:
                windows.append(win)
            return windows
        
        else:
            
            near_window = self.find_window(near,NEAR,ratio,span,result,thres)
            mid_window  = self.find_window(mid,MID,ratio,span,result,thres)
            far_window  = self.find_window(far,FAR,ratio,span,result,thres)
            
            for win in near_window[0]:
                windows.append(win)
            for win in mid_window[0]:
                windows.append(win)
            for win in far_window[0]:
                windows.append(win)
            
            thres_info = [[near_window[1],2000,10000],
                          [mid_window[1],666,2000],
                          [far_window[1],400,666]]
            
            return windows, thres_info    
        
        

    def find_window(self,nu,coef,ratio=0.3,span=100,result="nu",thres=False):
        
        window = []
        win = [0,0]
        
        max_signal = max(coef)
        min_signal = min(coef)
        
        signal_range = max_signal-min_signal
        
        threshold = min_signal + ratio * signal_range
        
        start = True
        
        for i,n in enumerate(nu):
        
                       
            if coef[i] < threshold and start == True:
                
                win[0] = n
                start = False
            
            if coef[i] > threshold and start == False:
                
                if n>=win[0]+span:
                    win[1] = n
                    start = True
                    window.append(np.array(win))
                else:
                    win[0] = n
                    start = True        
            """
            if i+1 == len(nu):
                if n>=win[0]+span:
                    win[1] = n
                    window.append(np.array(win))
            """      

        if result == "wav":
            wav_window = []
            for win in window[::-1]:
                low,high = "%.4g"%(10000./win[1]),"%.4g"%(10000./win[0])
                wav_window.append([float(low),float(high)])             
            window = wav_window

        if self.user_input["Save"]["Window"]["save"] == "true":
            with open(os.path.join(self.user_input["Save"]["Window"]["path"],
                                   self.user_input["Save"]["Window"]["name"]),"w") as f:
                for i in window:
                    f.write("%s-%s\n"%(i[0],i[1]))

        if thres:
            return window,threshold
        return window
     
    def spectra_window(self, nu, coef, type="A",threshold=200.,span=100.,min_signal=0,result="nu"):
        
        if type == "A":
        
            window = []
            
            start = True
            win = [0,0]
            
            self.stuff = []
            for i,n in enumerate(nu):
                
                if coef[i] < threshold and start == True:
                    win[0] = n
                    self.stuff.append(n)
                    start = False
                if coef[i] > threshold and start == False:
                    
                    self.stuff.append(n)
                    if n>=win[0]+span:
                        win[1] = n
                        start = True
                        window.append(np.array(win))
                    else:
                        win[0] = n
                        start = True
        
        elif type == "T":
            
            window = []
            start = True
            win = [0,0]
                        
            Min = min_signal
            Max = max(coef)#self.max_signal
            if threshold > 1:
                threshold = 1000/threshold
            
            
            threshold = Min+(Max-Min)*threshold
            self.threshold = threshold
            
            self.stuff = []
            for i,n in enumerate(nu):
                
                if coef[i] < threshold and start == True:
                    win[0] = n
                    self.stuff.append(n)
                    start = False
                if coef[i] > threshold and start == False:
                    
                    self.stuff.append(n)
                    if n>=win[0]+span:
                        win[1] = n
                        start = True
                        window.append(np.array(win))
                    else:
                        win[0] = n
                        start = True

        if self.user_input["Save"]["Window"]["save"] == "true":
            with open(os.path.join(self.user_input["Save"]["Window"]["path"],
                                   self.user_input["Save"]["Window"]["name"]),"w") as f:
                for i in window:
                    f.write("%s-%s\n"%(i[0],i[1]))
        
        
        
        wav_window = []
        
        for win in window[::-1]:
            low,high = "%.4g"%(10000./win[1]),"%.4g"%(10000./win[0])
            wav_window.append([low,high])
        
        
        if result == "nu":
            return window
        
        elif result == "wav":
            return wav_window
            
            
            
            

    def analyze_spectra_detection(self,nu,nu_window,trans,bio_trans,min_signal,method="max",result="bool"):
        """
        How to implement area under curve?
        bin size?
        """
        
        noise_level = 10
        comp = 2
        detection = False
        Detected = []
        Detected_num = []
        for i in nu_window:
            detected = False
            detected_num = 0
            reference =  trans[list(nu).index(i[0]):list(nu).index(i[1])]
            signal = bio_trans[list(nu).index(i[0]):list(nu).index(i[1])]
            
            
            # above certain ppm
            difference = max(signal-reference)*10**6
            # above certain comparision
            comparison = max((signal-min_signal)/(reference-min_signal))
        
            if difference > 3*noise_level:
                detection = True
                detected = True
                detected_num = 1
            if comparison > comp:
                detection = True
                detected = True
                detected_num = 1
                
            
            Detected.append(detected)
            Detected_num.append(detected_num)
        if result == "bool":
            return detection, Detected
        elif result == "num":
            return detection, Detected_num
        
 