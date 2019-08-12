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
Data Plotting Tools

More Modes is going to come

"""

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from matplotlib import ticker
ml = MultipleLocator(10)

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils as utils
from SEAS_Utils.common_utils.data_saver import check_path_exist, check_file_exist


class Plotter():
    """
    A more genralized plotter for multiple purpose needs
    """
    
    def __init__(self):
        pass
    
    
    def plot_xy(self, x, y):
        
        plt.plot(x,y)
        plt.show()


class Simulation_Plotter():
    """
    Plotter for the transit simulation that requires an user input cfg
    """
    '''
    def __init__(self, x, y,
                 line_style = "-", color = None,
                 title = "A plot", xlabel = "xlabel", ylabel = "ylabel",
                 show = True, time = None,
                 save = False, save_path = "", save_name = "Temp_Plot.png",
                 overwrite = False
                 ):

        """
        x and y are double arrays here... or in another plotter?
        or ... simply make this powerful?
        """
        
        self.Figure = plt.Figure()
        
        self.x = y
        self.y = y
        
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        
        self.line_style = line_style
        self.color = color
        
        

        self.show = show
        self.time = time
        
        
        self.save = save
        
        
        check_path_exist(save_path)
        
        self.savename = os.path.join(save_path,save_name)        
    '''
    
    def __init__(self, user_input):
        
        self.user_input = user_input
        
        fig_w = utils.to_int(self.user_input["Plotting"]["Figure"]["figsize_w"])
        fig_h = utils.to_int(self.user_input["Plotting"]["Figure"]["figsize_h"])
        
        self.fig = plt.figure(figsize=(fig_w, fig_h))
        
        self.ax = plt.gca()
        
        self.ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))

        
        xticks = ticker.MaxNLocator(10)
        self.ax.xaxis.set_major_locator(xticks) 
        
        
        self.x_unit  = self.user_input["Plotting"]["Figure"]["x_unit"]
        self.x_scale = self.user_input["Plotting"]["Figure"]["x_scale"]
        self.x_label = self.user_input["Plotting"]["Figure"]["x_label"]
        self.x_multi = utils.to_float(self.user_input["Plotting"]["Figure"]["x_multiplier"])
        
        self.y_unit  = self.user_input["Plotting"]["Figure"]["y_unit"]
        self.y_scale = self.user_input["Plotting"]["Figure"]["y_scale"]
        self.y_label = self.user_input["Plotting"]["Figure"]["y_label"]
        self.y_multi = utils.to_float(self.user_input["Plotting"]["Figure"]["y_multiplier"])
        
        self.Title   = self.user_input["Plotting"]["Figure"]["Title"]
        
        self.set_scale()
        plt.tick_params(axis='x', which='minor')
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        plt.title(self.Title)
        
    def plot_dot(self,x,y):
        
        plt.plot(x,y*self.y_multi,"*")
    
    
    def plot_hline(self,Y,Xlim,Label="HLine",Color = ""):
        
        if Color == "":
            
            plt.plot((10000./Xlim[1], 10000./Xlim[0]), (Y*self.y_multi, Y*self.y_multi), 'k-')
            
            #plt.axhline(y=Y*self.y_multi,xmin=,xmax=, linewidth=2, label=Label)
        else:
            plt.plot((10000./Xlim[1], 10000./Xlim[0]), (Y*self.y_multi, Y*self.y_multi), 'k-')
            #plt.axhline(y=Y*self.y_multi,xmin=10000./Xlim[1],xmax=10000./Xlim[0], linewidth=2, label=Label, color = Color)     
        
    
    def plot_xy(self, x, y, Label="Data", Color = "", Dtype="wn", Marker="-"):
        """
        assuming input is in wn... otherwise is um
        """
        if Dtype == "wn":
            if self.x_unit == "um":
                x = 10000./np.array(x)
            elif self.x_unit == "nm":
                x = 10000000./np.array(x)
            else:
                x = np.array(x)
        elif Dtype == "um":
            if self.x_unit == "wn":
                x = 10000./np.array(x)
            elif self.x_unit == "nm":
                x = np.array(x)*1000
            else:
                x = np.array(x)            
        
        y = np.array(y)
    
        if Color == "":
            if Marker == "-":
                plt_ref, = plt.plot(x*self.x_multi,y*self.y_multi,label=Label)   
            else:
                plt_ref, = plt.plot(x*self.x_multi,y*self.y_multi,label=Label, marker=Marker)            
        else:
            if Marker == "-":
                plt_ref, = plt.plot(x*self.x_multi,y*self.y_multi,label=Label,color=Color)   
            else:
                plt_ref, = plt.plot(x*self.x_multi,y*self.y_multi,label=Label,color=Color, marker=Marker)

        return plt_ref
  
    def plot_xy_list(self, xlist, ylist, ):
        pass
    
    
    def plot_bond_location(self,bond,max_signal,Color="r",Alpha=0.2):
        
        
        
        for k in bond:
            label = k[2]
            
            if self.x_unit == "um":
                if k[0] > k[1]:
                    up,down = 10000./k[0],10000./k[1]
                else:
                    up,down = 10000./k[1],10000./k[0]
                    
            elif self.x_unit == "nm":
                if k[0] > k[1]:
                    up,down = 10000000./k[0],10000000./k[1]    
                else:
                    up,down = 10000000./k[1],10000000./k[0]        
            elif self.x_unit == "wn":
                if k[0] > k[1]:
                    up,down = k[1],k[0]  
                else:
                    up,down = k[0],k[1]
            else:
                up,down = k[0],k[1]        
            
            plt.axvspan(up,down,facecolor=Color,alpha=Alpha)
            
            if k[2] == "C=C":
                max_signal -= 20*10**-6
                
            
            plt.text((up+down)/2, max_signal*self.y_multi, label,horizontalalignment='center')
    
    def plot_window(self, window, Color="k", Alpha=0.2):
        
        for k in window:
            if self.x_unit == "um":
                up,down = 10000./k[1],10000./k[0]
            elif self.x_unit == "nm":
                up,down = 10000000./k[1],10000000./k[0]          
            else:
                up,down = k[0],k[1]
            
            plt.axvspan(up,down,facecolor=Color,alpha=Alpha)

    def plot_bin(self,bin_centers, error_I, error_bar):
        
        
        print error_bar[0]*self.y_multi
        #sys.exit()
        
        plt.errorbar(bin_centers, error_I*self.y_multi, xerr=0, yerr=error_bar*self.y_multi, fmt='o',linewidth=2.0)

    def plot_SNR(self, bin_center, SNR):
        
        plt.plot(bin_center, SNR, ".")
    
    
    
    
        
    def plot_line(self):
        pass
    
    def set_scale(self):
        
        if self.x_scale == "log":
            self.ax.set_xscale('log')
        if self.y_scale == "log":
            self.ax.set_yscale('log')
    
    def set_legend(self, legends):
        plt.legend(handles=legends)

    def save_plot(self):
        
        save_dir  = self.user_input["Save"]["Plot"]["path"]
        save_name = self.user_input["Save"]["Plot"]["name"] 
        
        check_path_exist(save_dir)
        
        plt.savefig(os.path.join(save_dir,save_name))


    def save_window(self, nu_window):
        
        save_dir  = self.user_input["Save"]["Window"]["path"]
        save_name = self.user_input["Save"]["Window"]["name"]         

        check_path_exist(save_dir)
        with open(os.path.join(save_dir,save_name),"w") as file:
            for window in nu_window:
                file.write("%s-%s\n"%(window[0],window[1]))
            file.close()


    def show_plot(self):
        
        plt.show()  
    
