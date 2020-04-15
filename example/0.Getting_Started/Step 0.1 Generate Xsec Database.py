"""

To utilize the full capability of SEAS, 
the user will need to generate their own cross section database.

Run this code to generate. 


rewriting the old cross section generation methods to directly output into hdf5
Also clean up all old code with relation to cross section generation


Citation for Hapi:
R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016) [Link to article].



hdf5+gzip can compress cross section data as compared to raw npy outputs.
can further compress the data if less significant digits are used.
But how much is enough? certainly we don't need 10 digits but is 3 or 4 good enough?


Note that cross sections is usually exponentially increasing/decreasing
Take the log10 of the data further shrinks the datasize 

log10 + 4 digit significantly shrinks the datasize.

This approximation will yield a ~1% error on the cross section,
which is small compare to the errorbar on the data and errorbar of the abundance retrieved ... but for JWST?

further shrink using gzip compression with compressor_opt set to 9


--- What resolution /step should be chosen?
the difference between resolution is significant and higher compare to rounding error

but once it's binned down they're the same. So there's no benefit for going ultra-high resolution.
Just enough for what you're doing.

For JWST ~ 1000 dp, the xsec has ~10000 is fine enough.
which is why the wn_bin is different for each wavelength 
wn_bin = [[200,400,0.1],[400,2000,0.4],[2000,10000,2],[10000,30000,10]]

since JWST follows by wavelength and xsec is initially by wavenumber

hitran nu is different from np.arange(numin,numax,step). 
For standarization, normalize to np.arange().
This doesn't matter much for low res using JWST.



"""
import os
import sys
import h5py
from tqdm import trange
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.Common_Utils.configurable as config
import SEAS_Main.Cross_Section.Cross_Section_Calculator as csc

from SEAS_Main.Cross_Section.HITRAN_Match import HITRAN_Match


def test_generate_cross_section():
    """
    test cross section generation with simplistic inputs
    """
    # this is better, will be used for new xsec database
    wn_bin = [[200,400,0.1],[400,2000,0.4],[2000,10000,2],[10000,30000,10]]
    
    # this is to be consistent with old database
    wn_bin = [[400,2000,0.4],[2000,10000,2],[10000,30000,5]]
    
    
    nu_ = np.concatenate([np.arange(x[0],x[1],x[2]) for x in wn_bin])
    
    #T_Grid = csc.calculate_temperature_layers(T_Min=100, T_Max=800, Step=25)
    #P_Grid = csc.calculate_pressure_layers(P_surface = 1e5,P_Cutoff = 1e-5)
    user_input = config.Configuration("../../config/user_input.cfg")
    T_Grid = user_input["Xsec"]["Molecule"]["T_Grid"]
    P_Grid = user_input["Xsec"]["Molecule"]["P_Grid"]
    
    molecule = "HNO3"
    component = [HITRAN_Match[molecule],1,1]
    
    d_path = "../../SEAS_Input/Line_List/HITRAN_Line_List/%s"%molecule
    r_path = "../../SEAS_Input/Cross_Section/HDF5_New"


    P = 0.001
    T = 300
    for i in wn_bin:
        try:
            numin,numax,step = i
            xsec_calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step)
            nu,sigma = xsec_calc.hapi_calculator(P,T)
            plt.plot(nu,sigma)
            os.remove(os.path.join(d_path,"%s.header"%molecule))
            os.remove(os.path.join(d_path,"%s.data"%molecule))
            
            plt.show()
        except Exception as e:
            print(e)
        except:
            print("unknown")
            
def compress_cross_section():
    
    input = h5py.File("CO2.hdf5", "r")
    xsec_grid = input["results"]
    
    hdf5_store = h5py.File("CO2_new.hdf5", "w")
    hdf5_store.create_dataset("results", data=xsec_grid,compression="gzip",compression_opts=9)
    hdf5_store.close()

def examine_old_nu():

    DB_DIR = "../../SEAS_Input/Cross_Section/HDF5_DB"
        
    nu = h5py.File("%s/%s.hdf5"%(DB_DIR,"nu"), "r")
    nu = np.array(nu["results"])
    
    plt.plot(nu)
    
    wn_bin = [[200,400,0.1],[400,2000,0.4],[2000,10000,2],[10000,30000,10]]
    
    wn_bin = [[400,2000,0.4],[2000,10000,2],[10000,30000,5]]
    
    nu_ = np.concatenate([np.arange(x[0],x[1],x[2]) for x in wn_bin])
    
    plt.plot(nu_)
    plt.show()
        
def generate_cross_section_molecule(molecule):
    """
    generate cross section for a given molecule
    
    Known issue: using higher resolution hitran nu is not the same as np.arange(numin,numax,step).    
    """
    
    component = [HITRAN_Match[molecule],1,1]

    # this is better, will be used for new xsec database
    wn_bin = [[200,400,0.1],[400,2000,0.4],[2000,10000,2],[10000,30000,10]]
    
    # this is to be consistent with old database
    wn_bin = [[400,2000,0.4],[2000,10000,2],[10000,30000,5]]
    
    
    nu_ = np.concatenate([np.arange(x[0],x[1],x[2]) for x in wn_bin])
    

    #T_Grid = csc.calculate_temperature_layers(T_Min=100, T_Max=800, Step=25)
    #P_Grid = csc.calculate_pressure_layers(P_surface = 1e5,P_Cutoff = 1e-5)
    
    user_input = config.Configuration("../../config/user_input.cfg")
    T_Grid = np.array(user_input["Xsec"]["Molecule"]["T_Grid"],dtype=float)
    P_Grid = np.array(user_input["Xsec"]["Molecule"]["P_Grid"],dtype=float)/101300
    
    d_path = "../../SEAS_Input/Line_List/HITRAN_Line_List/%s"%molecule
    r_path = "../../SEAS_Input/Cross_Section/HDF5_New"

    # np.zeros maybe more adequate?
    xsec_grid = [[[] for y in range(len(P_Grid))] for x in range(len(T_Grid))]
    #sys.stdout = open('redirect.txt', 'w')
    for i in trange(len(wn_bin), desc='wn bin', leave=True):
        numin,numax,step = wn_bin[i]
        nu_std = np.arange(numin,numax,step)
        
        
        data = [[[] for y in range(len(P_Grid))] for x in range(len(T_Grid))]
        
        try:
            xsec_calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step)
            for i in trange(len(T_Grid),desc="Temperature"):
                T = T_Grid[i]
                for j in trange(len(P_Grid),desc="Pressure   "):
                    P = P_Grid[j]
                    #print("Generated %s xsec from %s to %s at T: %s and P: %s"%(molecule,numin,numax,P,T))
                    nu,sigma = xsec_calc.hapi_calculator(P,T) # note that P for hapi_calculator is in atm unit.
                    data[i][j] = sigma
                    
            
            
        except Exception as e:
            print("Missing HITRAN Data for %s from %s to %s"%(molecule,numin,numax))
            print(e)
            for i in trange(len(T_Grid),desc="Temperature"):
                T = T_Grid[i]
                for j in trange(len(P_Grid),desc="Pressure   "):
                    P = P_Grid[j]
                    data[i][j] = np.zeros(len(nu_std))
        except:
            print("unknown error")
            sys.exit()
        
        try:
            os.remove(os.path.join(d_path,"%s.header"%molecule))
            os.remove(os.path.join(d_path,"%s.data"%molecule))
        except:
            pass
        

        xsec_grid = np.concatenate([xsec_grid,data],2)
        
            
        #[float("%.4g"%x) for x in np.log10(data["results"][j][i])]
    
    
    r_path = ""
    
    hdf5_store = h5py.File(os.path.join(r_path,"%s.hdf5"%molecule), "w")
    hdf5_store.create_dataset("results", data=xsec_grid,compression="gzip",compression_opts=9)
    hdf5_store.close()
    
def generate_cross_section_database():
    """
    generate the entire cross section database that is used for atmosphere simulation
    """
    
    return


if __name__ == "__main__":
    
    #test_generate_cross_section()
    generate_cross_section_molecule("HNO3")
    #examine_old_nu()
    
    
    
    
    
    

