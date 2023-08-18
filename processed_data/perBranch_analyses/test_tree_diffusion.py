#! /usr/bin/env python

from treeswift import *
from math import sqrt

with open("../groundtruth_trees.txt","r") as fin:
    for line in fin:
        ID,treestr = line.strip().split()
        T = read_tree_newick(treestr)
        for node in T.traverse_preorder():
            if node.is_root():
                continue
            pnode = node.parent
            if pnode.label == "virtual":
                continue   
            x,y = [float(z) for z in node.node_params.split("(")[1][:-1].split(",")]
            x_par,y_par = [float(z) for z in pnode.node_params.split("(")[1][:-1].split(",")]
            t = node.edge_length
            print(ID,t,x,y,x_par,y_par,(x-x_par)/sqrt(t),(y-y_par)/sqrt(t))
