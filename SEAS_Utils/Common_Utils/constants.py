"""
Astrophysical and Atmospherical constants used in the simulation

all paramters, unless specified, are assumed to be in SI units


"""
import numpy as np
from scipy.constants import h,k,c




# Physical Parameters
HPlanck = 6.626070040e-34
c = 2.998e8
CLight = 2.998e8
V_c = 2.998e8
BoltK = 1.38064852e-23
mH = 1.6605389e-27
amu = 1.6605389e-27
G = 6.67e-11




R_Sun = 6.96e8
M_Sun = 2.0e30
T_Sun = 5770.0
L_Sun = 3.8e26

AU = 1.495e11
Psec = 206265*AU



# Planets

R_Earth = 6.4e6
M_Earth = 5.972e24
D_Earth = 1.5e11




M_Mercury = 3.285e23




# atmos
cm1_joules = 5.03445e22
AvoR = 6.0221415e23



HITRAN_Match = {"H2O":1,
                "CO2":2,
                "O3":3,
                "N2O":4,
                "CO":5,
                "CH4":6,
                "O2":7,
                "NO":8,
                "SO2":9,
                "NO2":10,
                "NH3":11,
                "HNO3":12,
                "OH":13,
                "N2":22,
                "H2O2":25,
                "C2H2":26,
                "C2H6":27,
                "HO2":33,
                "C2H4":38,
                "H2":45
                }

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
deg2rad=np.pi/180.
bar2barye=1.e6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
metricton2kg=1000. #1 metric ton in kg




