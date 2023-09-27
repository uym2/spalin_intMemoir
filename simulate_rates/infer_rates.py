#! /usr/bin/env python

from treeswift import *
import cvxpy as cp
import numpy as np
from sys import argv
from math import *

def label_to_frame(node_label):
    if node_label == "virtual":
        return 1
    frame,cellID = node_label.split("_")
    return int(frame)

input_file = argv[1]
smooth_factor = float(argv[2])
rate_factors = []
with open(argv[3],'r') as fin:
    for line in fin:
        rate_factors.append(float(line.strip()))

tree_list = read_tree_newick(input_file)
print(len(tree_list),len(rate_factors))

nframes = 216
nrates = nframes-1
A = []
b = []
for tree,r in zip(tree_list,rate_factors):
    for v in tree.traverse_preorder():
        v_frame = label_to_frame(v.label)
        if v.is_root():
            u_frame = 1
        else:
            u_frame = label_to_frame(v.parent.label)
        a = [0]*nrates        
        for i in range(u_frame-1,v_frame-1):
            a[i] = r
        A.append(a)
        b.append(v.edge_length*1000)  

for i in range(nrates-1):
    a = [0]*nrates
    a[i+1] = smooth_factor
    a[i] = -smooth_factor
    A.append(a)
    b.append(0)

A = np.array(A)
b = np.array(b)
x = cp.Variable(nrates,nonneg=True)
cost = cp.sum_squares(A @ x - b)
objective = cp.Minimize(cost)
prob = cp.Problem(objective)
prob.solve(solver=cp.GUROBI,verbose=False)
for i,y in enumerate(x.value):
    print(i,y)
