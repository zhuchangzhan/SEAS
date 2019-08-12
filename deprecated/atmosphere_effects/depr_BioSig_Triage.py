"""
This is a test to interpolate the NIST Spectra with Simulation format



change C to X
change spectra to something more general that covers Xsec, absorbance, Transmittance
move the interpolation to another module
"""

#from imports import *

import numpy as np
from scipy import interpolate
import os
import sys
import matplotlib.pyplot as plt

BoltK = 1.38*10**-23


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.jdx_Reader as jdx
import SEAS_Utils.common_utils.db_management2 as dbm

'''
def HITRAN_CIA():
    """
    steps of 25K from 200-3000
    each files contain 9981 data point. resolution is 1 cm^-1
    
    N2 contains 554... they are not the same. 
    """
    
    filename = "../../../input/CIA/H2-H2_2011.cia"
    #print util.check_file_exist(filename)
    
    nu,coef = [],[]
    temp = []
    t_nu, t_coef = [],[]
    start = True
    with open(filename,"r") as f:
        result = f.read().split("\n")
        for i in result:
            if i != "":
                try:
                    a,b = i.split()
                    t_nu.append(float(a))
                    t_coef.append(float(b))
                except:
                    if start:
                        start = False
                        temp.append(i.split()[4])
                    else:
                        temp.append(i.split()[4])
                        nu.append(t_nu)
                        coef.append(t_coef)
                        t_nu,t_coef = [],[]
            
    return nu,coef,temp    


def Exomol_spectra():
    filename = "../../../data/Exomol_xsec/PH3_300_test.txt"
    #print util.check_file_exist(filename)
    
    nu,coef = [],[]
    with open(filename,"r") as f:
        result = f.read().split("\n")
        for i in result:
            if i != "":
                a,b = i.split()
                nu.append(float(a))
                coef.append(float(b))
            
    return nu,coef
'''

def HITRAN_xsec():
    
    filename = "../../input/absorption_data/HITRAN_Cross_Section/O3/O3_300.0_0.0_29164.0-40798.0_04.xsc"
    nu,coef = [],[]
    with open(filename,"r") as f:
        result = f.read().split("\n")
        for i in result[1:]:
            if i != "":
                for j in i.split():
                    coef.append(j)
        
        numin = 29164.
        numax = 40798.
        
        wavmin = 10000./numax
        wavmax = 10000./numin
        
        npts = 5818
        
        wav = np.linspace(wavmin,wavmax,npts)
        
        nu = 10000./wav[::-1]
        
    return nu,coef


"""
def NIST_spectra(molecule,return_param):

    kwargs = {"db_name":"molecule_db.db",
          "user":"azariven",
          "dir":"/Users/mac/Workspace/BioSig/database",
          "DEBUG":False,"REMOVE":False,"BACKUP":False,"OVERWRITE":False}
    
    cross_db = dbm.database(**kwargs)   
    cross_db.access_db()   

    cmd = "SELECT inchikey From ID WHERE Formula='%s'"%molecule
    
    result = cross_db.c.execute(cmd)
    
    try:
        fetch = result.fetchall()[0][0]
    except:
        print "Molecule Doesn't Exist in NIST Database"
        sys.exit()
    
    path = os.path.join("/Users/mac/Workspace/BioSig/data/NIST_data/",fetch)
    
    filename = ""
    for j in os.listdir(path):
        if "jdx" in j:
            filename = os.path.join(path,j)
    
    data = jdx.JdxFile(filename)    
   
    
    if return_param[0] == "wl":
        x = data.wl()
    elif return_param[0] == "wn":
        x = data.wn()
        
    if return_param[1] == "T":
        y = data.trans()
    elif return_param[1] == "A":
        y = data.absorb()
    
    
    return x,y


"""
def HITRAN_spectra(molecule,spectra_param,return_param):
    
    
    kwargs = {"db_name":"cross_section_Sparse.db",
          "user":"azariven",
          "dir":"/Users/mac/Workspace/BioSig/database",
          "DEBUG":False,"REMOVE":False,"BACKUP":False,"OVERWRITE":False}
    
    cross_db = dbm.database(**kwargs)   
    cross_db.access_db()   
    
    unit = 0.0001
    
    Pref = 100000
    Tref = 300
    
    P = spectra_param[0]
    T = spectra_param[1]
    numin = spectra_param[2]
    numax = spectra_param[3]
    pathl = spectra_param[4]
    
        
    result = cross_db.c.execute("SELECT nu, coef FROM {} WHERE P={} AND T={} AND nu>={} AND nu<{} ORDER BY nu".format(molecule,P,T,numin,numax))
    fetch = np.array(result.fetchall()).T
    
    nu, coef = fetch
    n = Pref/(BoltK*Tref)
    absorb = n*coef*pathl*unit
    trans = np.exp(-absorb)
    
    
    if return_param[0] == "wl":
        x = 10000/nu[::-1]
        if return_param[1] == "T":
            y = trans[::-1]
        elif return_param[1] == "A":
            y = absorb[::-1]
        elif return_param[1] == "C":
            y = coef[::-1]
    elif return_param[0] == "wn":
        x = nu
        if return_param[1] == "T":
            y = trans
        elif return_param[1] == "A":
            y = absorb
        elif return_param[1] == "C":
            y = coef    
    return x,y



def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
        return idx-1#array[idx-1]
    else:
        return idx#array[idx]

def test_interpolate(x1,y1,x2, type):

    x1min = min(x1)
    x1max = max(x1)
    x2min = min(x2)
    x2max = max(x2)


    f = interpolate.interp1d(x1, y1)

    try:
        if x1min > x2min and x1max < x2max:
            print "A"
            left = find_nearest(x2,min(x1))+1
            right = find_nearest(x2,max(x1))
        
            if type == "A" or type == "C":
                yinterp_left = np.zeros(left)
                yinterp_right = np.zeros(len(x2)-right)
            elif type == "T":
                yinterp_left = np.ones(left)
                yinterp_right = np.ones(len(x2)-right)
            yinterp_middle = f(x2[left:right])
            yinterp = np.concatenate([yinterp_left,yinterp_middle, yinterp_right])
    
        elif x1min <= x2min and x1max < x2max:
            print "B"
            right = find_nearest(x2,max(x1))
            
            if type == "A" or type == "C":
                yinterp_right = np.zeros(len(x2)-right)
            elif type == "T":
                yinterp_right = np.ones(len(x2)-right)
            yinterp_middle = f(x2[:right])
            yinterp = np.concatenate([yinterp_middle, yinterp_right])
        
        elif x1min > x2min and x1max >= x2max:
            print "C"
            left = find_nearest(x2,min(x1))+1
        
            if type == "A" or type == "C":
                yinterp_left = np.zeros(left)
            elif type == "T":
                yinterp_left = np.ones(left)
            yinterp_middle = f(x2[left:])
            
            yinterp = np.concatenate([yinterp_left,yinterp_middle])
        
        else:
            print "D"
            yinterp = f(x2)
        
    except:
        print x1min,x1max
        print x2min,x2max
        sys.exit()
    
    
    return yinterp





def NIST_All_Spectra(x2):

    kwargs = {"db_name":"molecule_db.db",
          "user":"azariven",
          "dir":"/Users/mac/Workspace/BioSig/database",
          "DEBUG":False,"REMOVE":False,"BACKUP":False,"OVERWRITE":False}
    
    cross_db = dbm.database(**kwargs)   
    cross_db.access_db()   

    cmd = 'SELECT ID.smiles, ID.inchikey FROM ID,Spectra \
            WHERE ID.inchikey=Spectra.inchikey AND Spectra.has_spectra="Y"'
    
    result = cross_db.c.execute(cmd)
    data = np.array(result.fetchall()).T
    
    smiles = data[0]
    spectras = data[1]
    
    stuff = []
    for i,spectra in enumerate(spectras):
        path = os.path.join("/Users/mac/Workspace/BioSig/data/NIST_data/",spectra)

        filename = ""
        for j in os.listdir(path):
            if "jdx" in j:
                filename = os.path.join(path,j)
                break
        
        
        data = jdx.JdxFile(filename) 
        
        
        x = data.wn()
        y = data.absorb()
        
        yinterp = test_interpolate(x,y,x2,"Absorb")
        
        stuff.append(yinterp)
        
        """
        plt.plot(x,y)
        
        if i>10:
            break
        """

    """
    ax = plt.gca()
    ax.set_xscale('log')
    plt.tick_params(axis='x', which='minor')
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))  
    
    plt.show()
    """
    
    print len(stuff)
    return smiles, stuff

'''
def temperature_scaling(sigma, T_0, T):
    
    return sigma*np.sqrt(T_0/T)

def spectra_differences():
    P = 10000
    T = 300.
    numin = 400
    numax = 30000
    pathl = 0.002
    
    for pathl in []:
        molecule = "C2H2"
        bins = 1000
        x1,y1 = NIST_spectra(molecule,["wn","T"])
    
        x2,y2 = HITRAN_spectra(molecule,[P,T,numin,numax,pathl],["wn","T"])
    
        yinterp = test_interpolate(x1,y1,x2,"T")
    
    
        dp = len(x2)
        slices = np.linspace(0, dp, bins+1, True).astype(np.int)
        counts = np.diff(slices)
        mean_x_1 = np.add.reduceat(x2, slices[:-1]) / counts
        mean_y_1 = np.add.reduceat(yinterp, slices[:-1]) / counts
        
            
        
        dp = len(x2)
        slices = np.linspace(0, dp, bins+1, True).astype(np.int)
        counts = np.diff(slices)
        mean_y_2 = np.add.reduceat(y2, slices[:-1]) / counts
        
        
        plt.plot(mean_x_1,mean_y_2,"r")
        plt.plot(mean_x_1,mean_y_1,"b")
            
        plt.plot(mean_x_1,mean_y_2-mean_y_1,"g")
        
        print pathl, sum(mean_y_2-mean_y_1)/len(mean_y_2)

'''

if __name__ == "__main__":
    
    """
    Pref = 10000.
    Tref = 300.
    nref = Pref/(BoltK*Tref)
    
    lref = 0.05

    unit = 10000

    numin = 400
    numax = 30000
    pathl = 1000000

    molecule = "H2"    
    
    x1,y1 = HITRAN_xsec()
    x2,y2 = HITRAN_spectra(molecule,[100000,300,numin,numax,pathl],["wn","T"])
    
    yinterp = test_interpolate(x1,y1,x2,"C")
    
    plt.plot(10000./x1,y1)
    plt.plot(10000./x2,yinterp)
    plt.show()
    
    sys.exit()
    """
    
    
    """
    Pref = 10000.
    Tref = 300.
    nref = Pref/(BoltK*Tref)
    
    lref = 0.05

    unit = 10000

    numin = 400
    numax = 30000
    pathl = 1000000

    molecule = "H2"
    
    
    
    
    a,b = HITRAN_spectra(molecule,[100000,300,numin,numax,pathl],["wl","T"])
    plt.plot(a,b)
    plt.show()
    
    sys.exit()
    
    
    x1,y1 = NIST_spectra(molecule,["wn","T"])
    x2,y2 = HITRAN_spectra(molecule,[100000,300,numin,numax,pathl],["wn","T"])

    
    y1 = np.array(y1)+(1-(np.mean(y1)+np.median(y1))/2)
    y1new = []
    for i in y1:
        if i > 1:
            y1new.append(1)
        else:
            y1new.append(i)
    y1 = y1new
    
    yinterp = test_interpolate(x1,y1,x2,"T")
    sigma = -np.log(yinterp)/(nref*lref)*unit
    """
    
    
    
    """
    for P in [100000,10000,1000,100]:
        for T in [250,300,350]:
            n = P/(BoltK*T)
            
            y = np.exp(-n*sigma*lref*0.0001)
            plt.plot(x2,y)
    
    """
    """
    
    
    x2,sigma,temp = HITRAN_CIA()
    
    #plt.plot(x1,y1)
    #plt.plot(10000./x2,yinterp)
    for i in range(10):
        
        plt.plot(x2[0],sigma[i],label=temp[i])
    plt.title("H2 CIA")
    plt.xlabel("wavenumber")
    plt.ylabel("intensity")
    plt.legend()
    plt.show()
    
    
    """
    
    
    
    
    

    """
    x1,y1 = HITRAN_spectra(molecule,[100000,300,numin,numax,pathl],["wl","T"])
    
    x2,y2 = HITRAN_spectra(molecule,[10000,300,numin,numax,pathl],["wl","T"])
    
    x3,y3 = HITRAN_spectra(molecule,[1000,300,numin,numax,pathl],["wl","T"])

    #yinterp = test_interpolate(x1,y1,x2,"T")

    plt.plot(x1,y1,"r")
    plt.plot(x2,y2,"g")
    plt.plot(x3,y3,"b")

    plt.show()



    """
    """
    x1,y1 = NIST_spectra("C2H2")
    x2,y2 = HITRAN_spectra("CO")
    #NIST_All_Spectra(x2)

    """







