import numpy as np



def multi_column_file_loader(path,spliter=None,type="float",skip=0):
    """
    load data from files that contains multiple columns
    """
    
    with open(path) as data_file:   
        data = [x.split(spliter) for x in data_file.read().split("\n")[skip:]]
        if data[-1] == [] or data[-1] == None or data[-1] == "":
            data = data[:-1]
        
        
        data = [list(x) for x  in zip(*data)]    
        
        if type == "float":
            return np.array(data,dtype=np.float)
        elif type == "int":
            return np.array(data,dtype=np.int)
        elif type == "mixed":
            return data
