#! /usr/bin/env python

from treeswift import *
from sys import argv

tree = read_tree_newick(argv[1])
spatial_file = argv[2]
S = {}
with open(spatial_file,'r') as fin:
    for line in fin:
        name,x,y = line.strip().split(",")
        S[name] = (x,y)

for node in tree.traverse_preorder():
    if not node.is_root():
        v = node.label
        x_v,y_v = S[v]
        h_v = v.split("_")[0]
        u = node.parent.label
        x_u,y_u = S[u]
        h_u = u.split("_")[0]        
        print(x_u,y_u,h_u,x_v,y_v,h_v)
