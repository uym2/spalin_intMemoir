#! /bin/bash

for x in s*c*; do 
    echo $x
    python ~/my_gits/problin/run_problin.py -c $x/characters.txt -t $x/without_location.nwk -o $x/without_location --delimiter comma -p prior.csv -L "0 0"
    head -n1 $x/without_location_annotations.txt | nw_distance -mp -sa -n - | sort > $x/without_location_nMutations.txt 
done    
