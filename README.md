# SEAS
Simulator for Exoplanet Atmosphere Spectra (SEAS)
-- Exoplanet Atmosphere Benchmark/validation/sanity check tool

# Installation
Code download and installation using command line:

1. git clone https://github.com/zhuchangzhan/SEAS.git
2. virtualenv SEAS_Env
3. Source SEAS_Env/Bin/Activate
4. pip install requirements.txt

# Testing the code:
All tested execution files are under /bin
Some simple functionalities and code demos will work without cross section database.
However, to use the full functionality of the code, including generating your own spectra,
will require creating your own absorption cross section database. 
The code for this is under Example/...

# Molecular Absorption Cross Section Database
Absorption cross section are pre-generated and stored in SEAS_Input/ as hdf5 files


# Boxcar Transmission Model
example/Atmosphere_Simulation_Template/Example_Boxcar_TS_Model.py

The boxcar transmission model is a overly simplified transmission spectra model which calculates the transmission of light passing though a 1m standard temperature and pressure box that contains a defined set of molecules.

The purpose of this model is to construct a toy example for ease of implementing and testing the retrieval model using MCMC or Machine Learning.

Retrieval on the boxcar transmission model will try to retrieve the moleculear composition of the boxcar. Understanding the retrieval of this model will provide useful insights for implementing retrieval on the actual atmosphere model.