from treeswift import *
from random import gauss
from math import *
import numpy as np
from scipy.optimize import minimize,LinearConstraint 
from random import shuffle

def sep_gauss_locations(tree,leaf_location):
# assume Gaussian diffusion and opposite separation forces
    idx = 0
    for node in tree.traverse_preorder():
        if not node.is_leaf():
            node.idx = idx
            idx += 1
            c1,c2 = node.children
            c1.sign = 1
            c2.sign = -1            
    for node in tree.traverse_leaves():
        node.idx = idx
        idx += 1        
    N = idx + len(list(tree.traverse_leaves()))-1
    s_base_idx = idx
    A_x = []
    b_x = []
    A_y = []
    b_y = []
    for node in tree.traverse_postorder():
        if node.is_root():
            a_xy = [0]*N
            a_s = [0]*N
            c1,c2 = node.children
            a_xy[node.idx] = 1/c1.edge_length + 1/c2.edge_length
            a_xy[c1.idx] = -1/c1.edge_length
            a_xy[c2.idx] = -1/c2.edge_length
            a_xy[s_base_idx+node.idx] = c1.sign/c1.edge_length + c2.sign/c2.edge_length
            a_s[node.idx] = 1/c1.edge_length + 1/c2.edge_length
            a_s[c1.idx] = -1/c1.edge_length
            a_s[c2.idx] = -1/c2.edge_length
            a_s[s_base_idx+node.idx] = 1/c1.edge_length + 1/c2.edge_length + 0.5
            A_x.append(a_xy)
            A_x.append(a_s)
            A_y.append(a_xy)
            A_y.append(a_s)
            b_x.append(0)
            b_x.append(0)
            b_y.append(0)
            b_y.append(0)        
        elif node.is_leaf():
            a = [0]*N
            a[node.idx] = 1
            A_x.append(a)
            A_y.append(a)
            b_x.append(leaf_location[node.label][0])
            b_y.append(leaf_location[node.label][1])
        else:
            p = node.parent
            c1,c2 = node.children
            a_xy = [0]*N
            a_s = [0]*N
            a_xy[node.idx] = 1/node.edge_length + 1/c1.edge_length + 1/c2.edge_length
            a_xy[c1.idx] = -1/c1.edge_length
            a_xy[c2.idx] = -1/c2.edge_length
            a_xy[p.idx] = -1/node.edge_length
            a_xy[s_base_idx+node.idx] = c1.sign/c1.edge_length + c2.sign/c2.edge_length
            a_xy[s_base_idx+p.idx] = -node.sign/node.edge_length
            a_s[node.idx] = 1/c1.edge_length + 1/c2.edge_length
            a_s[c1.idx] = -1/c1.edge_length
            a_s[c2.idx] = -1/c2.edge_length
            a_s[s_base_idx+node.idx] = 1/c1.edge_length + 1/c2.edge_length + 0.5
            A_x.append(a_xy)
            A_x.append(a_s)
            A_y.append(a_xy)
            A_y.append(a_s)
            b_x.append(0)
            b_x.append(0)
            b_y.append(0)
            b_y.append(0)        
    x = np.dot(np.linalg.pinv(A_x),b_x)
    y = np.dot(np.linalg.pinv(A_y),b_y)
    #print(x[-10:])
    for node in tree.traverse_postorder():
        node.location = (x[node.idx],y[node.idx])
    return x,y

def gaussian_locations(tree,leaf_location):
# assume Gaussian diffusion
    idx = 0
    for node in tree.traverse_preorder():
        node.idx = idx
        idx += 1
    N = idx
    A_x = []
    b_x = []
    A_y = []
    b_y = []
    for node in tree.traverse_postorder():
        if node.is_root():
            a = [0]*N
            #a[node.idx] = 1
            #a[node.children[0].idx] = -1
            c1,c2 = node.children
            a[node.idx] = 1/c1.edge_length + 1/c2.edge_length
            a[c1.idx] = -1/c1.edge_length
            a[c2.idx] = -1/c2.edge_length
            A_x.append(a)
            A_y.append(a)
            b_x.append(0)
            b_y.append(0)
        elif node.is_leaf():
            a = [0]*N
            a[node.idx] = 1
            A_x.append(a)
            A_y.append(a)
            b_x.append(leaf_location[node.label][0])
            b_y.append(leaf_location[node.label][1])
        else:
            p = node.parent
            c1,c2 = node.children
            a = [0]*N
            a[node.idx] = 1/node.edge_length + 1/c1.edge_length + 1/c2.edge_length
            a[c1.idx] = -1/c1.edge_length
            a[c2.idx] = -1/c2.edge_length
            a[p.idx] = -1/node.edge_length
            A_x.append(a)
            A_y.append(a)
            b_x.append(0)
            b_y.append(0)        
    x = np.dot(np.linalg.pinv(A_x),b_x)
    y = np.dot(np.linalg.pinv(A_y),b_y)
    for node in tree.traverse_postorder():
        node.location = (x[node.idx],y[node.idx])
    return x,y

def average_locations(tree,leaf_location):
    # estimate the location of each node simply by a weighted average of the children locations
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.location = leaf_location[node.label]
        else:
            C = node.children
            C_x = [c.location[0] for c in C]    
            C_y = [c.location[1] for c in C]    
            C_w = [1/c.edge_length for c in C]    
            x = sum(a*w for a,w in zip(C_x,C_w))/sum(C_w)
            y = sum(a*w for a,w in zip(C_y,C_w))/sum(C_w)
            node.location = (x,y)

def t_locations(tree,leaf_location,df=13,scale=2.5):
    def nllh(x,*args): # the negative log-likelihood function (without constants)
        coordinate = args[0] # either 0 or 1
        nllh = 0
        for node in tree.traverse_preorder():
            if node.parent is not None:
                xc = leaf_location[node.label][coordinate] if node.is_leaf() else x[node.idx]
                xp = x[node.parent.idx]
                t = node.edge_length
                dx = (xc-xp)/sqrt(t)/scale
                nllh += (df+1)/2*log(1+dx**2/df)
        return nllh        
    gaussian_locations(tree,leaf_location) # the gaussian estimated locations are decorated on the tree using attribute "location"
    # idex internal nodes, ignore leaf nodes
    idx = 0
    x0 = []
    y0 = []
    for node in tree.traverse_preorder():
        if node.is_leaf():
            node.idx = None
        else:    
            node.idx = idx
            x0.append(node.location[0])
            y0.append(node.location[1])
            idx += 1
    results_x = minimize(nllh,x0,args=(0),method="CG")
    x = results_x.x
    results_y = minimize(nllh,y0,args=(1),method="CG")
    y = results_y.x
    for node in tree.traverse_postorder():
        node.location = (x[node.idx],y[node.idx])
    return x,y

def main():
    with open("../groundtruth_trees.txt","r") as fin:
        for line in fin:
            ID,treestr = line.strip().split()
            # the estimation: estimate using the true leaf locations
            T = read_tree_newick(treestr)
            leaf_location = {}
            for node in T.traverse_leaves():    
                x,y = [float(z) for z in node.node_params.split("(")[1][:-1].split(",")]
                leaf_location[node.label] = (x,y) 
            #gaussian_locations(T,leaf_location)
            #average_locations(T,leaf_location)
            sep_gauss_locations(T,leaf_location)
            d2root = {}
            b2root = {}
            for node in T.traverse_preorder():
                node.b2root = node.parent.b2root + 1 if not node.is_root() else 0-int(node.label == "virtual")
                node.d2root = node.parent.d2root + node.edge_length if not node.is_root() else 0
                d2root[node.label] = node.d2root
                b2root[node.label] = node.b2root
                if not node.is_leaf() and node.label != "virtual":
                    true_x,true_y = [float(z) for z in node.node_params.split("(")[1][:-1].split(",")]
                    est_x,est_y = node.location
                    print(ID,node.label,d2root[node.label],b2root[node.label],"estimated",true_x,true_y,est_x,est_y,sqrt((est_x-true_x)**2+(est_y-true_y)**2))
            # the permutation test: permutate the leaf labels to get the baseline
            '''nrep = 10
            labels_perm = [node.label for node in T.traverse_leaves()]
            for r in range(nrep):
                T1 = read_tree_newick(treestr)
                # shuffle
                shuffle(labels_perm)
                i = 0
                for node in T1.traverse_leaves():
                    node.label = labels_perm[i]
                    i += 1
                # estimate the location of the shuffled tree    
                gaussian_locations(T1,leaf_location)
                for node in T1.traverse_preorder():
                    if not node.is_leaf() and node.label != "virtual":
                        true_x,true_y = [float(z) for z in node.node_params.split("(")[1][:-1].split(",")]
                        est_x,est_y = node.location
                        print(ID,node.label,d2root[node.label],b2root[node.label],"permutated",true_x,true_y,est_x,est_y,sqrt((est_x-true_x)**2+(est_y-true_y)**2))'''
if __name__ == "__main__":
    main()        
