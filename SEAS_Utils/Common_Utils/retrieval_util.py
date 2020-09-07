"""

Collection of useful functions and utils when conducting retrieval analysis 
and useful notes when doing retrieval using HMC and MCMC









"""
import os,sys,time
import numpy as np


import pymc3 as pm
import theano.tensor as tt



def sample(model, cache_dir = "",cache_file="temp_trace_file.trace",resample=False, **kwargs):
    """
    wrapper of the pymc3 sample to allow loading trace from save.
    """
    
    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)
    
    filepath = os.path.join(cache_dir,cache_file)
    
    if resample or not os.path.exists(filepath):
        trace = pm.sample(**kwargs)
        pm.save_trace(trace, directory=filepath,overwrite=True)
    else:
        trace = pm.load_trace(filepath, model=model)
        
    return trace

def default_trace():
    """
    trace = pm.sample
    
    tune=2000, 
    draws=2000, 
    chains=2, 
    cores=2, 
    step=xo.get_dense_nuts_step()
        
    """
    
def useful_exoplanet_functions():
    """
    Useful exoplanet functions

    get_dense_nuts_step()
    
    # useful to first optimize the parameters to find the “maximum a posteriori” (MAP) parameters 
    # and then start the sampler from there. 
    # This is useful here because MCMC is not designed to find the maximum of the posterior; 
    # it’s just meant to sample the shape of the posterior.
    # map is also the mode of the posterior distribution
    map_soln = optimize() 
    
    
    # Better than spamming Deterministic nodes ?
    with model:
        ...
        mod = func(x,*args)
        eval_in_model(mod,map_soln)
    
    with model:
        mod2 = func(xgrid,*args)
        eval_in_model(mod2,map_soln)
    
    with model:
        y_grid = a0 + a1 * x_grid + a2 * x_grid ** 2
        for i, sample in enumerate(get_samples_from_trace(trace, size=50)):
            samples[i] = eval_in_model(y_grid, sample)
        
    get_samples_from_trace(trace, size=50)
    """
    pass

def useful_pymc3_functions():
    """
    
    """
    
    pass

def pymc3_distributions():

    """
    pm.Normal()
    
    # PyMC3 includes the construct Bound for placing constraints on existing probability distributions. 
    # It modifies a given distribution to take values only within a specified interval.
    pm.Bound(<distribution>,lower=, upper=)
    
    
    """
    
    
    pass

def pymc3_extra():
    """
    # Add the likelihood
    pm.Normal("obs", mu=mod, sd=pm.math.exp(logs), observed=y)
    
    most of the time you don't really care about obs, it's just to link the observation with model?
    
    
    once the model is created, we can evaluate the model using xo.eval_in_model to inspect the model
    and whether it's well created or not.
    This can be done before sampling happens
    
    """
    return

