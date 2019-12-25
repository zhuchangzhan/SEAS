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
Download the HITRAN Collision Induced Absorption Database

The HITRAN CIA exist as txt files.
use web scrapper to download the data

As of July 2017, the X-X_2016.cia files are not found on the website and can't be scraped. 
Consequently we will download the X-X_2011.cia datafiles for now.



"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.data_downloader.web_scraper as ws
from SEAS_Utils.common_utils.DIRs import HITRAN_CIA


if __name__ == "__main__":
    
    scraper = ws.HITRAN_CIA_crawler("http://hitran.org/cia/",HITRAN_CIA)
    scraper.download()



