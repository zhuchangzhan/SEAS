"""

This is the interpolation function that needs optimization
The runtime for grid_interpolate takes about 1second, 
ideally I want to get it down to 10ms, if that's even possible.

We can potentially store the data better, have less decimals, etc 
but the bulk of the time is in RegularGridInterpolator.

Perhaps we can write our own numba function for parallel calculation
since it's being called 12000 times? (each wavelength is interpolated though the 24x9 grid)

These are just ideas, feel free to take it away how you like it
"""


import sys
import time
import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator as rgi


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms'%(method.__name__, (te - ts) * 1000))
        return result
    return timed

@timeit
def grid_interpolate():
    
    DB_DIR = "data"
    
    # Precalculated Grids
    P_Grid = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]
    T_Grid = [100,150,200,250,275,300,325,350,400]
    
    # Cross section data for CH4. 
    # wavelength is the "x" for plot, 
    # xsec_grid  is the "y" for plot
    # xsec_grid  is a (24,9,12000) array.  
    wavelength = np.array(h5py.File("%s/%s.hdf5"%(DB_DIR,"nu"), "r")["results"])
    xsec_grid  = np.array(h5py.File("%s/%s.hdf5"%(DB_DIR,"CH4"), "r")["results"])
    
    # Actual measurements of pressure and temperature we want to interpolate to
    normalized_pressure    = [91483.44,76196.15,62819.36,51203.54,41203.06,32676.27,25519.61,19650.63,14928.96,11241.97,8449.276,6350.326,4772.792,3587.15,2696.038,2026.294,1522.927,1144.605,860.2655,646.5602,485.9432,365.2267,274.4979,206.3077,155.0572,116.5383,87.58818,65.82972,49.47645,37.18566,27.94808,21.00528,15.7872,11.86538,8.917816,6.702472,5.037459,3.786069,2.845542,2.138659,1.607378,1.208077,0.90797,0.6824141,0.5128903,0.3854792,0.2897193,0.2177482,0.1636556,0.1230007,0.09244523,0.06948019,0.05222008,0.03924768]
    normalized_temperature = [280.786111,266.358333,251.930556,237.502778,223.075,208.647222,196.405694,186.593333,178.404306,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0]
    
    normalized_xsec = np.zeros((54,12000))
    fn = rgi((np.log10(P_Grid[::-1]),T_Grid, wavelength),xsec_grid[::-1])
    for i,(P_E,T_E) in enumerate(zip(np.log10(normalized_pressure),normalized_temperature)):
        pts = np.array([np.ones(len(wavelength))*P_E,
                        np.ones(len(wavelength))*T_E,
                        wavelength]).T
        normalized_xsec[i] = fn(pts) 
            
    return normalized_xsec


if __name__ == "__main__":
    xsec = grid_interpolate()