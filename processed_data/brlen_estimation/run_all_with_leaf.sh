#! /bin/bash

for x in s*c*; do 
    echo $x
    python run_spalin_leaf_locations.py $x > $x/spalin_output_leaf.log 2>&1
done    
