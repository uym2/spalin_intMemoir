#!/bin/bash

# Preprocess the character matrices and groundtruth trees
# 1. rename the cells in the matrices to match that of the trees
# 2. prune the tree to the same leaf set as the character matrices
# 3. put each pair of (tree,matrix) to a folder named by the experiment ID

data_source="/Users/uym2/Documents/intMEMMOIR_Data/Figure3/reconstructionParameters/filteredData"
temp=`mktemp`

# create subfolders and put in the data matrices
for x in $data_source/*data.txt; do
    y=`basename $x | sed -e "s/_data.txt//g" -e "s/_//g"`
    mkdir $y
    cat $x | awk '{print $1;}' | sed -e "s/^/216_/g" -e "s/216_cell/cell/g" > $temp
    cat $x | awk '{print $2;}' | sed -e "s/0/,0/g" -e "s/1/,1/g" -e "s/2/,2/g" -e "s/state/,state/g" | paste $temp - | awk '{print $1$2;}' > $y/characters.txt
done    

rm $temp

# get the trees and locations
while read line; do
    ID=`echo $line | awk '{print $1;}'`
    echo $ID
    python extract_tree.py $line `tail -n+2 $ID/characters.txt | awk -F "," '{print $1;}'`
    grep "^216" $ID/true_locations.txt > $ID/leaf_locations.txt
done < ../groundtruth_trees.txt
