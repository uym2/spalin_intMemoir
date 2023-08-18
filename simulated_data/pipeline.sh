#! /bin/bash

#bash copy_data.sh
#for x in s*; do
#    echo $x
#    python ~/my_gits/problin/simulate_sequence.py -k 10 --alphabetsize 2 -t $x/true_tree.nwk -r 5 --mu 0.006 -p $x/k10
#done    
for r in {1..5}; do
    for x in s*; do
        echo $x $r
        echo spalin leaves
        python ~/my_gits/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk -S $x/leaf_locations.txt --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_spalin_leaves.txt --nInitials 10 > $x/k10_r$r\_spalin_leaves.log 2>&1
        echo spalin all
        python ~/my_gits/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk -S $x/true_locations.txt --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_spalin_all.txt --nInitials 10 > $x/k10_r$r\_spalin_all.log 2>&1
        echo problin
        python ~/my_gits/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_problin.txt --nInitials 10 > $x/k10_r$r\_problin.log 2>&1
    done 
done    
