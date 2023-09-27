#! /bin/bash

#bash copy_data.sh

#for x in s[0-9]*c[0-9]*; do
#    echo $x
#    python scale_tree.py $x
#    python ~/my_gits/problin/simulate_sequence.py -k 10 --alphabetsize 2 -t $x/true_scaled_tree_const.nwk -r 1 --mu 1 -p $x/k10_const
#    python ~/my_gits/problin/simulate_sequence.py -k 10 --alphabetsize 2 -t $x/true_scaled_tree_lin.nwk -r 1 --mu 1 -p $x/k10_lin
#    python ~/my_gits/problin/simulate_sequence.py -k 10 --alphabetsize 2 -t $x/true_scaled_tree_exp.nwk -r 1 --mu 1 -p $x/k10_exp

#    python count_zeros.py $x/k10_const_character_matrix.csv
#    python count_zeros.py $x/k10_lin_character_matrix.csv
#    python count_zeros.py $x/k10_exp_character_matrix.csv
#done   

for x in s[0-9]*c[0-9]*; do
    echo $x
    #python scale_tree_exp2.py $x
    python ~/my_gits/problin/simulate_sequence.py -k 200 --alphabetsize 2 -t $x/true_scaled_tree_exp2.nwk -r 1 --mu 1 -p $x/k200_exp2
    python ~/my_gits/problin/simulate_sequence.py -k 10 --alphabetsize 2 -t $x/true_scaled_tree_exp2.nwk -r 1 --mu 1 -p $x/k10_exp2
done    

#for x in s[0-9]*c[0-9]*; do
#    echo $x
#    python ~/my_gits/problin/simulate_sequence.py -k 200 --alphabetsize 2 -t $x/true_scaled_tree_const.nwk -r 1 --mu 1 -p $x/k200_const
#    python ~/my_gits/problin/simulate_sequence.py -k 200 --alphabetsize 2 -t $x/true_scaled_tree_lin.nwk -r 1 --mu 1 -p $x/k200_lin
#    python ~/my_gits/problin/simulate_sequence.py -k 200 --alphabetsize 2 -t $x/true_scaled_tree_exp.nwk -r 1 --mu 1 -p $x/k200_exp

#    python count_zeros.py $x/k200_const_character_matrix.csv
#    python count_zeros.py $x/k200_lin_character_matrix.csv
#    python count_zeros.py $x/k200_exp_character_matrix.csv
#done   

