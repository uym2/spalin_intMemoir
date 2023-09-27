#! /bin/bash

for x in s*; do
    for type in const lin exp exp2; do
        echo $x $type
        python ~/my_gits/problin/run_problin.py -c $x/k10_$type\_character_matrix.csv -t $x/true_tree.nwk --delimiter comma -p $x/k10_$type\_priors.csv -o $x/k10_$type\_problin --nInitials 10 --ultrametric > $x/k10_$type\_problin.log 2>&1
    done
done 
