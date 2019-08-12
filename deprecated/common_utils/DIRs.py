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
This file contains some of the common directories used in the code

Maybe this should be a config file as well? The User should be allow to edit this file.

Config file can also reduce the amount of os.path.join needed.... a simple conversion?

This is temporary, will rework in 0.7 to be more robust
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))
ROOT = sys.path[0]


Temp_DIR = os.path.join(ROOT,"output/Temp_Result")

Data_DIR = os.path.join(ROOT,"input")
Molecule_Absorption = os.path.join(Data_DIR,"absorption_data")
HITRAN_CIA          = os.path.join(Molecule_Absorption,"HITRAN_CIA")
NIST_Spectra        = os.path.join(Molecule_Absorption,"NIST_Spectra")
Exomol_Xsec         = os.path.join(Molecule_Absorption,"Exomol_Cross_Section")
Stellar_Spectra     = os.path.join(Molecule_Absorption,"Stellar_Spectra")
HITRAN_Lines        = os.path.join(Molecule_Absorption,"HITRAN_Line_List")
HITRAN_Water_Lines  = os.path.join(HITRAN_Lines, "H2O")
HITRAN_PH3_Lines    = os.path.join(HITRAN_Lines, "PH3")
HITRAN_CH4_Lines    = os.path.join(HITRAN_Lines, "CH4")
HITRAN_CH3OH_Lines    = os.path.join(HITRAN_Lines, "CH3OH")
HITRAN_Molecule_List= os.path.join(HITRAN_Lines, "HITRAN_Molecule_List.txt")

Aerosol_Data        = os.path.join(Data_DIR,"aerosol_data")
HITRAN_Mineral      = os.path.join(Aerosol_Data,"HITRAN_mineral")

Atmosphere_Data     = os.path.join(Data_DIR,"atmosphere_data")
TP_Profile_Data     = os.path.join(Atmosphere_Data, "TP_Profile")
Mixing_Ratio_Data   = os.path.join(Atmosphere_Data, "Mixing_Ratio")

DB_DIR = os.path.join(Data_DIR,"database")

Example_DB    = os.path.join(DB_DIR,"Example")
Demo_DB       = os.path.join(DB_DIR,"Demo")
Simulation_DB = os.path.join(DB_DIR,"Simulation_Band")
Exomol_DB     = os.path.join(DB_DIR,"Exomol")

molecule_info = os.path.join(Data_DIR, "molecule_info")

Intermediate_DIR = os.path.join(Data_DIR,"temporary")

Result_DIR = Data_DIR = os.path.join(ROOT,"result")


