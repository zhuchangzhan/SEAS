"""

# will be used for new xsec database
wn_bin = [[200,400,0.1],[400,2000,0.4],[2000,10000,2],[10000,30000,10]]

d_path = "../../SEAS_Input/Line_List/HITRAN_Line_List/%s"%molecule
r_path = "../../SEAS_Input/Cross_Section/HDF5_New"

There is an error with no linelist data (HNO3) that is not addressed in this demo
Fixes is implemented in cross section database generation.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.Cross_Section.Cross_Section_Calculator as csc
from SEAS_Main.Cross_Section.HITRAN_Match import HITRAN_Match


def test_calculate_cross_section(molecule="CO2",
                                 d_path = "",r_path = "",
                                 P=1,T=300.,
                                 numin=400.,numax=2000.,step=0.4,
                                 SL_Flag=False):
    """
    calculate cross section with specific inputs
    
    Parameters
    ----------
    molecule : str
        Name of the molecule. Has to be one which HITRAN recognize
    d_path : str
        filepath for where the linelist data will be stored
    r_path : str
        filepath for where the calculated cross section will be stored
    P : float
        Pressure [bar]
    T : float
        Temperature [K]
    numin : float
        minimum wavenumber [cm^-1]
    numax : float
        maximum wavenumber [cm^-1]
    step : float
        wavenumber increment step size. 
    SL_Flag : bool
        Flag for whether the downloaded hitran linelist will be saved or not. 
    
    
    Returns
    -------
    """
    
    try:
        component = [HITRAN_Match[molecule],1,1]
    except:
        print("molecule: %s not recognized, not in HITRAN?"%molecule)


    xsec_calc = csc.cross_section_calculator(d_path,molecule,component,numin,numax,step)
    nu, sigma = xsec_calc.hapi_calculator(P,T)
    
    plt.plot(nu,sigma)
    
    if not SL_Flag:
        os.remove(os.path.join(d_path,"%s.header"%molecule))
        os.remove(os.path.join(d_path,"%s.data"%molecule))
        
    
    plt.show()

def test_calculate_cross_section_multi(molecule="CO2",
                                       d_path = "",r_path = "",
                                       P=1,T=300.,
                                       wn_bin=[[400,2000,0.4],[2000,10000,2],[10000,30000,5]],
                                       SL_Flag=False):
    """
    Calculate cross section with specific inputs for multiple wavenumber bins.
    This function is useful because small step size for shorter wavelength is not necessary.
    Most parameters similar to function above with exception to wrapping numin, numax, and step in a nested list.
    
    Parameters
    ----------
    wn_bin : list of lists
        list of wavenumber bins. Each sublist contains [numin, numax, step]
        
    
    """
    
    try:
        component = [HITRAN_Match[molecule],1,1]
    except:
        print("molecule: %s not recognized, not in HITRAN?"%molecule)

    for i in wn_bin:
        numin,numax,step = i
        xsec_calc = csc.cross_section_calculator(d_path,molecule,component,numin,numax,step)
        nu, sigma = xsec_calc.hapi_calculator(P,T)
    
        plt.plot(nu,sigma)
    
    if not SL_Flag:
        os.remove(os.path.join(d_path,"%s.header"%molecule))
        os.remove(os.path.join(d_path,"%s.data"%molecule))
        
    
    plt.show()


if __name__ == "__main__":
    
    
    test_calculate_cross_section(SL_Flag=True)
    
    
    
    
    