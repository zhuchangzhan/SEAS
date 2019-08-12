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
Web downloader



"""
import os
import urllib
from SEAS_Utils import isfile

def get_HITRAN_ID():
    
    file = open("../../input/database/Molecule_Info/HITRAN_Molecule_List.txt").read().split("\n")
    ID = []
    for i,molecule in enumerate(file):
        ID.append([molecule,i+1,1])
    return ID


class downloader():
    
    def __init__(self,url,path):
        self.url = url
        self.path = path
        self.download()
    
    def download(self):
        urllib.urlretrieve(self.url,os.path.join(self.path,self.url.split("/")[-1]))
        

class HITRAN_Line_List_downloader():
    """
    Note that as of July 2017, the following molecules can not be downloaded from HITRAN
    and must be handled manually by going to 
    https://www.cfa.harvard.edu/hitran/HITRAN2012/HITRAN2012/Supplemental/
    
    which you will find .par files to the following molecules
    
    30 SF6 Sulfur Hexafluoride
    35 ClONO2 Chlorine Nitrate
    42 CF4 Carbon Tetrafluoride
    """
    
    def __init__(self,
                 outpath = "../../input/database/absorption_data/User_Defined",
                 molecule = ["H2O",1,1],
                 numin = 0,
                 numax = 50000,
                 down = True,
                 VERBOSE = True):
    
        self.n = molecule[0]
        self.m = molecule[1]
        self.i = molecule[2]
        self.min = numin
        self.max = numax
        
        self.path = "/".join([outpath,self.n])
        
        if down:
            filedir = "".join([self.path,"/",self.n,".data"])
            if not isfile(filedir):
                self.download()
            else:
                print self.n,
                print " data already exist, no new data downloaded"
                if VERBOSE:
                    print "Note that the current checker is not aware of data completeness"
                    print "If you're downloading data for previously not covered"
    
    def download(self):
        
        import SEAS_Aux.cross_section.hapi as hp
        
        hp.db_begin(self.path)
        try:
            hp.fetch(self.n,self.m,self.i,self.min,self.max)
        except KeyError:
            print self.n, " load error, data not downloaded"
        except:
            print self.n, " load error, data not downloaded"
    

class HITRAN_CIA_downloader():
    
    def __init__(self):
        pass
    
    
class NIST_Spectra_downloader():
    
    
    def __init__(self):








    