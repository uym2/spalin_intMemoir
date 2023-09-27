#! /usr/bin/env python

from problin_libs.sequence_lib import read_sequences
from problin_libs.SpaLin_solver import SpaLin_solver
from sys import argv
from treeswift import *
import os

k = 10
m = 3
Q = []

for i in range(k):
    q = {j+1:1/m for j in range(m)}
    q[0] = 0
    Q.append(q)

data = argv[1] 
print("dataset",data)
msa,_ = read_sequences(data+"/characters.txt",filetype="charMtrx")
tree_obj = read_tree_newick(data+"/true_tree.nwk")
T = tree_obj.newick()    

all_locations = {}
with open(data+"/leaf_locations.txt",'r') as fin:
#with open(data+"/true_locations.txt",'r') as fin:
    for line in fin:
        cellID,x,y = line.strip().split(",")
        all_locations[cellID] = (float(x),float(y))

D = {'charMatrices':[msa]}
prior = {'Q':[Q]}
params = {'nu':0,'phi':0}

sigma = 22.470982031241125 # pre-estimated from data
params['sigma'] = sigma
D['locations'] = all_locations 
#mySolver = SpaLin_solver(msa,Q,T,leaf_locations,sigma)
mySolver = SpaLin_solver([T],D,prior,params)
mySolver.optimize_locations()
for x in mySolver.inferred_locations:
    print(x,mySolver.inferred_locations[x])
