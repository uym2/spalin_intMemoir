#! /bin/bash

for x in s*c*; do 
    echo $x
    python run_spalin.py $x > $x/spalin_output.log 2>&1
done    
