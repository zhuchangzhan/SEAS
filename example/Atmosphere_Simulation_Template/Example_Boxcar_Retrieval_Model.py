"""
Here is a construction site for building SEAS.
Will have to clean up once finished
Let's build the user input section first
"""

import os
import sys
import tqdm
import hashlib
from scipy import stats
import matplotlib.pyplot as plt

#ROOT = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../..')
sys.path.insert(0, "../..")    

from SEAS_Utils.Common_Utils.constants import *
import SEAS_Utils.Common_Utils.configurable as config
import SEAS_Utils.System_Utils.optimization as opt
import SEAS_Utils.Common_Utils.data_loader as load 
import SEAS_Main.Physics.astrophysics as calc 
import SEAS_Main.Simulation.transmission_spectra_simulator as TS
from SEAS_Main.Physics.noise import Photon_Noise

import pymc3 as pm
import corner  # https://corner.readthedocs.io
import theano
import theano.tensor as tt
from theano.compile.ops import as_op

user_input = config.Configuration("../../config/user_input.cfg")
VERBOSE = bool(user_input["Data_IO"]["Logging"]["VERBOSE"])


#@as_op(itypes=[tt.lscalar], otypes=[tt.lscalar])
def generate_spectra(bin_edges,user_input,logH2O,logCH4):
    
    
    H2O = pm.math.exp(logH2O)
    CH4 = pm.math.exp(logCH4)

    # use percent for consistency for now
    Profile = {}
    Profile["Molecule"] = {}
    
    # even though it's only 1 layer, still use list to keep code consistent
    Profile["Temperature"]      = [300]
    Profile["Pressure"]         = [100000]
    Profile["Molecule"]["H2O"]  = [H2O]
    Profile["Molecule"]["CH4"]  = [CH4]
    Profile["Molecule"]["N2"]   = [100 - H2O - CH4]

    # Loading TP_Profile, MR_Profile, Molecule List, 
    user_input = load.load_Atmosphere_Profile(user_input, scenario_file=Profile)

    # Load Atmosphere model and generate theoretical spectra
    simulation = TS.Transmission_Spectra_Simulator(user_input)
    simulation.load_boxcar_model()
    user_input = simulation.user_input
    
    #print(user_input["Spectra"]["Total_Transit_Signal"])
    
    x       = user_input["Spectra"]["Wavelength"]
    true_y  = user_input["Spectra"]["Total_Transit_Signal"]
    
    print(type(true_y))
    
    sys.exit()
    
    
    
    true_y_bin, bin_edges, binnumber = stats.binned_statistic(10000./x[::-1], true_y[::-1], bins=bin_edges)

    
    return true_y_bin

def generate_spectra_std(bin_edges,user_input,logH2O,logCH4):
    
    
    H2O = np.exp(logH2O)
    CH4 = np.exp(logCH4)

    # use percent for consistency for now
    Profile = {}
    Profile["Molecule"] = {}
    
    # even though it's only 1 layer, still use list to keep code consistent
    Profile["Temperature"]      = [300]
    Profile["Pressure"]         = [100000]
    Profile["Molecule"]["H2O"]  = [H2O]
    Profile["Molecule"]["CH4"]  = [CH4]
    Profile["Molecule"]["N2"]   = [100 - H2O - CH4]

    # Loading TP_Profile, MR_Profile, Molecule List, 
    user_input = load.load_Atmosphere_Profile(user_input, scenario_file=Profile)

    # Load Atmosphere model and generate theoretical spectra
    simulation = TS.Transmission_Spectra_Simulator(user_input)
    simulation.load_boxcar_model()
    user_input = simulation.user_input
    
    x       = user_input["Spectra"]["Wavelength"]
    true_y  = user_input["Spectra"]["Total_Transit_Signal"]
    
    true_y_bin, bin_edges, binnumber = stats.binned_statistic(10000./x[::-1], true_y[::-1], bins=bin_edges)

    
    return true_y_bin


@opt.timeit
def Retrieval_Boxcar_Model_Architecture():
    
    global user_input
    
    user_input["Prototype"]["Source"]="Boxcar"
    
    # use percent for consistency for now
    Profile = {}
    Profile["Molecule"] = {}
    
    # even though it's only 1 layer, still use list to keep code consistent
    Profile["Temperature"]      = [300]
    Profile["Pressure"]         = [100000]
    Profile["Molecule"]["H2O"]  = [0.1]
    Profile["Molecule"]["CH4"]  = [0.1]
    Profile["Molecule"]["N2"]   = [100 - Profile["Molecule"]["H2O"][0] - Profile["Molecule"]["CH4"][0]]


    # Loading TP_Profile, MR_Profile, Molecule List, 
    user_input = load.load_Atmosphere_Profile(user_input, scenario_file=Profile)
    
    # Load absorption cross section for all molecule and effect (CIA, Cloud, etc)
    user_input = load.load_Absorption_Cross_Section(user_input,False)

    
    # Load Atmosphere model and generate theoretical spectra
    simulation = TS.Transmission_Spectra_Simulator(user_input)
    simulation.load_boxcar_model()
    user_input = simulation.user_input
    
    true_logs = np.log(0.01)
    
    x       = user_input["Spectra"]["Wavelength"]
    true_y  = user_input["Spectra"]["Total_Transit_Signal"]
 
    
    noise = Photon_Noise(user_input)
    bin_edges, bin_width, bin_centers = noise.determine_bin()  
    
    true_y_bin, bin_edges, binnumber = stats.binned_statistic(10000./x[::-1], true_y[::-1], bins=bin_edges)

    x_bin = bin_centers
    
    
    y_noise = true_y_bin + np.exp(true_logs)*np.random.randn(len(bin_centers))
    
    
    plt.plot(x_bin,y_noise,".k",markersize=1)
    plt.plot(10000./x,true_y,linewidth=1,alpha=0.5)
    
    
    plt.xscale("log")
    plt.show()
    
    sys.exit()
    
    
    with pm.Model() as model:
        
        logH2O = pm.Uniform("logH2O", lower=np.log(1e-8), upper=np.log(100))#, testval=np.log(1))
        logCH4 = pm.Uniform("logCH4", lower=np.log(1e-8), upper=np.log(100))#, testval=np.log(1))
        logs   = pm.Uniform("logs", lower=-5, upper=5)
        
        obs = pm.Normal("obs", 
                        mu=generate_spectra(bin_edges,user_input,logH2O,logCH4), 
                        sd=pm.math.exp(logs), 
                        observed=y_noise)

    
        # This is how you will sample the model. Take a look at the
        # docs to see that other parameters that are available.
        trace = pm.sample(draws=1000, tune=1000, chains=2, cores=2)
        
    
    pm.traceplot(trace, varnames=["logH2O", "logCH4", "logs"])
    plt.show()
    
    print(pm.summary(trace, varnames=["logH2O", "logCH4", "logs"]))
    
    samples = pm.trace_to_dataframe(trace, varnames=["logH2O", "logCH4", "logs"])
    corner.corner(samples, truths=[0.1, 0.1, true_logs]);
    plt.show()


    plt.figure(figsize=(7, 7))


    for i in np.random.randint(len(trace) * trace.nchains, size=50):
        
        xlogH2O,xlogCH4,xlogs = trace["logH2O"][i],trace["logCH4"][i],trace["logs"][i]
        
        plt.plot(x_bin, generate_spectra_std(bin_edges,user_input,xlogH2O,xlogCH4)+ np.exp(xlogs)*np.random.randn(len(x_bin)), color="C0", lw=1, alpha=0.3)
    
    plt.plot(x_bin, y_noise, 'x', label='data')
    plt.plot(x_bin, true_y_bin + np.exp(true_logs)*np.random.randn(len(x_bin)), label='true regression line', lw=1., c='r')
    plt.title('Posterior predictive regression lines')
    plt.legend(loc=0)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xscale("log")
    plt.show()



if __name__ == "__main__":
    
    Retrieval_Boxcar_Model_Architecture()
















