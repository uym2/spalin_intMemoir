#! /usr/bin/env python

from treeswift import *
from sys import argv

tree = read_tree_newick(argv[1])

for node in tree.traverse_preorder():
    b = node.edge_length if node.edge_length is not None else 0
    print(node.label,b)

