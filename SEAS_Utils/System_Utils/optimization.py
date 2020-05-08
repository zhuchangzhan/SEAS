import sys
import time
import h5py
import inspect

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

def name_of_var(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def print_this(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    print(str([k for k, v in callers_local_vars if v is var][0])+': '+str(var))

@timeit
def memory_check(variable):
    
    checker = """Memory: %s"""%sys.getsizeof(variable)
    
    #print(name_of_var(checker))
    #print_this(checker)
    return checker 

@timeit
def check_arg_kwarg(*arg,**kwarg):
    print(arg)
    print(kwarg)
    print(*arg)
    #print(**kwarg)

def compress_h5py_database(input_file,output_file,tablename="results",compress="gzip",compress_opts=9):
    """
    reduce the storage space of uncompressed database
    
    """
    
    input = h5py.File(input_file, "r")
    input_data = input[tablename]
    
    hdf5_store = h5py.File(output_file, "w")
    hdf5_store.create_dataset(tablename, 
                              data=input_data,
                              compression=compress,
                              compression_opts=compress_opts)
    hdf5_store.close()

if __name__ == "__main__":
    你大爷 = ["yes","no","1"]
    item = 1
    you = [1,2,3,4]
    
    dicto = {"1":1,"2":2}
    check_arg_kwarg(你大爷,item,you,dicto,a=1,b=2,c=3)
    
    
    
    