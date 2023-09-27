#! /bin/bash

r=1
#for r in {1..5}; do
    for x in s*; do
        echo $x $r
        echo problin
        python ~/my_gits/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_problin_ultrametric --nInitials 10 --ultrametric > $x/k10_r$r\_problin_ultrametric.log 2>&1
    done 
#done    
