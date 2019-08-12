"""

tools to manipulte data files

includes
    bin, trim, stitch, etc
    
    also include interpolation stuff
"""

import numpy as np
from scipy import stats
from scipy import interpolate

def trim_data(xlist,ylist,up,down):
    
    
    for i,info in enumerate(xlist):
        if info > up:
            start = i
            break
    for j,info in enumerate(xlist[i:]):
        if info > down:
            end = j+i
            break
    return xlist[start:end],ylist[start:end]


def bin_data(xlist,ylist,bin_num=100, bin_range=(0,1)):

    info = stats.binned_statistic(xlist, ylist,bins=bin_num, 
                                  statistic='mean',range=bin_range)
    
    bin_means, bin_edges, binnumber = info 
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2

    return bin_centers,bin_means


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
        return idx-1#array[idx-1]
    else:
        return idx#array[idx]

def interpolate_data(x1,y1,x2, Type):

    x1min = min(x1)
    x1max = max(x1)
    x2min = min(x2)
    x2max = max(x2)


    f = interpolate.interp1d(x1, y1)
    
    if x1min > x2min and x1max < x2max:
        #print "A"
        
        left = find_nearest(x2,min(x1))+1
        right = find_nearest(x2,max(x1))
    
        if Type == "A" or Type == "C":
            yinterp_left = np.zeros(left)
            yinterp_right = np.zeros(len(x2)-right)
        elif Type == "T":
            yinterp_left = np.ones(left)
            yinterp_right = np.ones(len(x2)-right)
        yinterp_middle = f(x2[left:right])
        yinterp = np.concatenate([yinterp_left,yinterp_middle, yinterp_right])

    elif x1min <= x2min and x1max < x2max:
        #print "B"
        right = find_nearest(x2,max(x1))
        
        if Type == "A" or Type == "C":
            yinterp_right = np.zeros(len(x2)-right)
        elif Type == "T":
            yinterp_right = np.ones(len(x2)-right)
        yinterp_middle = f(x2[:right])
        yinterp = np.concatenate([yinterp_middle, yinterp_right])
    
    elif x1min > x2min and x1max >= x2max:
        #print "C"
        left = find_nearest(x2,min(x1))+1
    
        if Type == "A" or Type == "C":
            yinterp_left = np.zeros(left)
        elif Type == "T":
            yinterp_left = np.ones(left)
        yinterp_middle = f(x2[left:])
        
        yinterp = np.concatenate([yinterp_left,yinterp_middle])
    
    else:
        #print "D"
        yinterp = f(x2)

    
    
    return yinterp

def merge_list(list_d, param=5):
    """
    merge similar x in an array
    
    """
    
    new_list = []
    compare_list = []
    tag = []
    prev = -100
    begin = True
    
    if list_d == []:
        return []
    elif len(list_d) == 1:
        return list_d
    else:
        pass
    
    for i in list_d:
        
        if begin:
            begin = False
            prev = i
            continue
    
        compare_list.append(abs(i-prev) <= param)
        prev = i
    
 
    
    for j,id in enumerate(compare_list):
    
        if id == True:
            if tag == []:
                tag.append(list_d[j])
                tag.append(list_d[j+1])
            else:
                tag.append(list_d[j+1])
        
        else:
            try:
                value = tag[0]+(tag[-1]-tag[0])/2
                new_list.append(value)
            except:
                new_list.append(list_d[j])
            tag = []
    
    if tag != []:
        new_list.append(tag[0]+(tag[-1]-tag[0])/2)
    
    if compare_list[-1] == False:
        new_list.append(list_d[-1])

        
    
    return new_list





def list_merger(indexs, value, param=5, method = "max"):
    
    new_list = []
    compare_list = []
    tag = []
    prev = -100
    begin = True
    
    if indexs == []:
        return []
    elif len(indexs) == 1:
        return indexs
    else:
        pass

     
    for i in indexs:
        
        if begin:
            begin = False
            prev = i
            continue
    
        compare_list.append(abs(i-prev) <= param)
        prev = i 
    
    
    for j,id in enumerate(compare_list):
        if id == True:
            if tag == []:
                tag.append(indexs[j])
                tag.append(indexs[j+1])
            else:
                tag.append(indexs[j+1])
        else:
            if tag == []:
                tag.append(indexs[j])
            
            
            if method == "max":
                y = map(value.__getitem__, tag)
                max_index = list(value).index(max(y))
                new_list.append(max_index)
            elif method == "mid":
                new_list.append(tag[0]+(tag[-1]-tag[0])/2)
            tag = []    
    
    if tag != []:
        if method == "max":
            y = map(value.__getitem__, tag)
            max_index = list(value).index(max(y))
            new_list.append(max_index)
        elif method == "mid":
            new_list.append(tag[0]+(tag[-1]-tag[0])/2)
        
        
    if compare_list[-1] == False:
        new_list.append(indexs[-1])

    
    return new_list  
    
