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
Web Crawler to download data from websites


"""

import os
import urllib2
import httplib
from bs4 import BeautifulSoup


from SEAS_Utils import to_float
from SEAS_Utils.common_utils.configurable import ConfigurableObject
import web_downloader as wd


# example crawler class
class web_crawler(ConfigurableObject):
    
    def __init__(self):
        
        pass

    def lookup(self):
        
        pass
    
    def save(self):
        
        pass


class nist_crawler(web_crawler):
    
    def __init__(self):
        pass
    


class HITRAN_CIA_crawler():
    
    def __init__(self,url,path):
        
        self.url = url
        self.path = path
    
    
    def download(self):

        page = urllib2.urlopen(self.url).read()
        soup = BeautifulSoup(page,"lxml")
        data = soup.findAll("li")
        
        for i in data:
            
            if "CIA" in i.a["href"]:
                link = "".join(["http://hitran.org",i.a["href"]])
                if "2016" in link:
                    link = link.replace("2016","2011")
                
                print "downloading %s"%link
                wd.downloader(link,self.path)
        
        print "All data downloaded to %s"%self.path
                
                
                
                
    
    
    
    