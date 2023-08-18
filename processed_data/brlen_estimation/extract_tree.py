#/usr/bin/env python

from treeswift import *
from sys import argv

ID = argv[1]
tree = read_tree_newick(argv[2])
S = argv[3:]

tree1 = tree.extract_tree_with(S,suppress_unifurcations=True)
tree1.write_tree_newick(ID + "/true_tree.nwk")

extracted_cells = set(node.label for node in tree1.traverse_preorder() if node.label != "virtual")
with open(ID + "/true_locations.txt",'w') as fout:
    for node in tree.traverse_preorder():
        if node.label in extracted_cells:
            location =  node.node_params.split("(")[1][:-1]
            fout.write(node.label + "," + location + "\n")

