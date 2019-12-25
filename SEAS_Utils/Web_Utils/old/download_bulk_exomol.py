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

Download cross section data from exomol website


Note that exomol cross sections doesn't have temperature < 296, also it doesn't have pressure broadening



"""

import os
import mechanize
import numpy as np
from bs4 import BeautifulSoup

from SEAS_Utils.common_utils.data_downloader import link_loader
from SEAS_Utils.common_utils.DIRs import Exomol_Xsec

def lookup(link,data):

    # initiate the mechanize and read the relevent link
    br = mechanize.Browser()
    br.set_handle_robots(False)
    br._factory.is_html = True
    br = link_loader(link,mecha=br,mechanize=True)
    
    # set parameters for cross section
    br.select_form(nr=1)
    br.form.set_value(str(data["id_dnu"]),id="id_dnu")
    br.form.set_value(str(data["id_numin"]),id="id_numin")    
    br.form.set_value(str(data["id_numax"]),id="id_numax")    
    br.form.set_value(str(data["id_T"]),id="id_T")       
    br.form.find_control(type="checkbox").items[0].selected = data["checkbox"]
    
    # read output link and find the cross section data link
    result = br.submit()
    page = result.read()
    soup = BeautifulSoup(page, "html.parser")
    
    
    links = [(link["href"] if ".sigma" in link["href"] else "") for link in soup.findAll("a")]
    for i in links:
        if i!="":
            datalink = "http://exomol.com"+i
            name = i.split("/")[-1]
            break
        
    # read cross section data from datalink    
    moresoup = BeautifulSoup(link_loader(datalink), "html.parser")
    
    return moresoup, name
    

        
    
    
        
        
if __name__ == "__main__":
    
    
    for T in [100,150,200,250,275,300,325,350,400]:
        
        if T< 296:
            Te = 296
        else:
            Te = T
        
        molecule = "31P-1H3"
        link = "http://exomol.com/xsec/%s/"%molecule
        
        data = {}
        data["id_dnu"]   = 1
        data["id_numin"] = 400
        data["id_numax"] = 9999
        data["id_T"]     = Te
        data["checkbox"] = True
        
        result,name = lookup(link,data)
        
        
        lines = str(result).split("\n")
        xdata, ydata = [],[]
        for line in lines:
            splitted = line.split()
            if splitted == []:
                pass
            else:
                xdata.append(splitted[0])
                ydata.append(splitted[1])
        
        if T<296:
            name = name.replace("296","%s"%T)
            
    
        d_path = os.path.join(Exomol_Xsec, molecule,name.replace("sigma","npy"))
    
        np.save(d_path,[xdata,ydata])
        print name


        
        
        
        
        