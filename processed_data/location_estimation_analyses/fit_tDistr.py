# /usr/bin/env python
from math import sqrt
import scipy.stats

with open("tree_diffusion.txt","r") as fin:
    X = []
    Y = []
    fin.readline() # skip header
    for line in fin:
        _,_,_,_,_,_,dx_norm,dy_norm = line.strip().split()
        X.append(float(dx_norm))
        Y.append(float(dy_norm))
        
print(scipy.stats.t.fit(X))
print(scipy.stats.t.fit(Y))
