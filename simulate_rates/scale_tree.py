#! /bin/bash

from sys import argv
from treeswift import *
from math import *

wdir = argv[1]
mu = 0.006

# constant rate
tree_const = read_tree_newick(wdir+"/true_tree.nwk")
for node in tree_const.traverse_preorder():
    if node.edge_length is None:
        node.edge_length = 0
    node.edge_length *= mu
tree_const.write_tree_newick(wdir+"/true_scaled_tree_const.nwk") 

# linearly degraded rate
mu_max = 0.008
mu_min = 0.004
nrates = 215
delta_mu = (mu_max-mu_min)/(nrates-1)
tree_lin = read_tree_newick(wdir+"/true_tree.nwk")

for node in tree_lin.traverse_preorder():
    if node.is_root():
        node.mu_start = mu_max
    else:
        node.mu_start = node.parent.mu_end
    mu = node.mu_start
    el = 0
    if node.edge_length is None:
        node.edge_length = 0
    for i in range(int(node.edge_length)):
        el += mu
        mu -= delta_mu
    node.edge_length = el
    node.mu_end = mu

tree_lin.write_tree_newick(wdir+"/true_scaled_tree_lin.nwk")

# exponentially degraded rate
mu_max = 0.008
mu_min = 0.004
nrates = 215
log_delta_mu = (log(mu_max)-log(mu_min))/(nrates-1)
tree_lin = read_tree_newick(wdir+"/true_tree.nwk")

for node in tree_lin.traverse_preorder():
    if node.is_root():
        node.mu_start = mu_max
    else:
        node.mu_start = node.parent.mu_end
    mu = node.mu_start
    el = 0
    if node.edge_length is None:
        node.edge_length = 0
    for i in range(int(node.edge_length)):
        el += mu
        mu = exp(log(mu)-log_delta_mu)
    node.edge_length = el
    node.mu_end = mu

tree_lin.write_tree_newick(wdir+"/true_scaled_tree_exp.nwk")
