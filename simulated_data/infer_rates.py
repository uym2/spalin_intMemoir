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

w = float(argv[1])
input_file = "k10_r1_problin_ultrametric_trees.nwk"
tree_list = read_tree_newick(input_file)

nframes = 216
nrates = nframes-1
A = []
b = []
for tree in tree_list:
    for v in tree.traverse_preorder():
        v_frame = label_to_frame(v.label)
        if v.is_root():
            u_frame = 1
        else:
            u_frame = label_to_frame(v.parent.label)
        a = [0]*nrates        
        for i in range(u_frame-1,v_frame-1):
            a[i] = 1
        A.append(a)
        b.append(v.edge_length)  

for i in range(nrates-1):
    a = [0]*nrates
    a[i+1] = w
    a[i] = -w
    A.append(a)
    b.append(0)

A = np.array(A)
b = np.array(b)
x = cp.Variable(nrates,nonneg=True)
cost = cp.sum_squares(A @ x - b)
objective = cp.Minimize(cost)
#constraints = [np.zeros(nrates)+0.004 <= x, x <= np.zeros(nrates)+0.008]
#prob = cp.Problem(objective,constraints)
prob = cp.Problem(objective)
prob.solve(solver=cp.MOSEK,verbose=False)
#print(sqrt(np.mean([(y-0.006)**2 for y in x.value])))
for i,y in enumerate(x.value):
    print(i,y)
