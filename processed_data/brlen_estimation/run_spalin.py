#! /usr/bin/env python

from problin_libs.sequence_lib import read_sequences
from problin_libs.ML_solver import ML_solver
from problin_libs.EM_solver import EM_solver
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
# remove branch lengths
#for node in tree_obj.traverse_preorder():
#    node.edge_length = None
T = tree_obj.newick()    

all_locations = {}
#with open(data+"/leaf_locations.txt",'r') as fin:
with open(data+"/true_locations.txt",'r') as fin:
    for line in fin:
        cellID,x,y = line.strip().split(",")
        all_locations[cellID] = (float(x),float(y))

D = {'charMtrx':msa}
prior = {'Q':Q}
params = {'nu':0,'phi':0}
'''print("Optimization without using locations")
mySolver = ML_solver([T],D,prior,params)
mySolver.optimize(initials=20,fixed_nu=1e-8,fixed_phi=1e-8,verbose=False,ultra_constr=True)
mySolver.trees[0].write_tree_newick(data+"/without_location.nwk")'''

print("Optimization using locations")
sigma = 22.470982031241125 # pre-estimated from data
params['sigma'] = sigma
D['locations'] = all_locations 
#mySolver = SpaLin_solver(msa,Q,T,leaf_locations,sigma)
mySolver = SpaLin_solver([T],D,prior,params)
mySolver.optimize(initials=1,fixed_nu=0,fixed_phi=0,optimize_brlens=False,verbose=True,ultra_constr=True)
#mySolver.trees[0].write_tree_newick(data+"/with_all_location.nwk")
mySolver.trees[0].write_tree_newick("temp")
'''with open(data+"/with_location.log",'w') as fout:    
    fout.write("sigma: " + str(mySolver.params.sigma) + "\n")
    for cell in mySolver.inferred_locations:
        x,y = mySolver.inferred_locations[cell]
        fout.write(cell + " " + str(x) + " " + str(y) + "\n") '''
