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

Data Downloader 

"""

import time
import urllib2
import httplib

def link_loader(link, safeguard=True, iter=3, mecha="", mechanize=False,timer=0):
    """
    load a url into page information
    if saveguard is true, will try to compensite for errors occuring during failed load
    """
    
    if timer!=0:
        time.sleep(timer)
    
    if iter == 0:
        print "iteration link loading failed"
        return ""
    
    if not safeguard:
        if mechanize:
            mecha.open(link)
            return mecha
        else:
            return urllib2.urlopen(link)
    
    try:
        if mechanize:
            mecha.open(link)
            return mecha
        else:
            return urllib2.urlopen(link)
        
    except httplib.BadStatusLine:
        print "BadStatusLine, attempting to resolve"
        time.sleep(10)
        response = link_loader(iter=iter-1, mechanize=mechanize)
            
    except httplib.IncompleteRead:
        print "IncompleteRead, attempting to resolve"
        time.sleep(10)
        response = link_loader(iter=iter-1, mechanize=mechanize)   
               
    except KeyboardInterrupt:
        return ""
    

def search():
    pass
    
    





