import numpy as np

def myfunc(a, b):
    print(a,b)
    "Return a-b if a>b, otherwise return a+b"
    if a > b:
        return [a - b]
    else:
        return [a + b]
    



vectfunc = np.vectorize(myfunc,otypes=[list],cache=False)
print(list(vectfunc([1,2,3],[1,2,3])))