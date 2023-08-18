#! /bin/bash

for x in s*c*; do
    sed -e "s/,/ /g" s10c1/true_locations.txt | sort -k1,1 > temp1
    grep -v "sigma" s10c1/with_leaf_location.log | sort -k1,1 > temp2
    join temp1 temp2 | grep -v "^216" | awk '{print $1,"brlen_coestimated",$2,$3,$4,$5,sqrt(($4-$2)**2+($5-$3)**2);}' |sed -e "s/^/$x /g"
    rm temp1 temp2
done    
