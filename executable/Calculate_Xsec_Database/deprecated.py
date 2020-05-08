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