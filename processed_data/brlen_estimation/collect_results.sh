#! /bin/bash

for x in s*/with_leaf_location.nwk; do echo $x; python collect_brlen.py $x | sort -k1,1 > `dirname $x`/`basename $x .nwk`_brlen.txt; done
for x in s*/with*brlen.txt; do y=`echo $x | sed -e "s/\// /g" -e "s/_brlen.*//g"`; join `dirname $x`/true_tree_brlen.txt $x | join `dirname $x`/true_tree_nodeAge.txt - | sed -e "s/^/$y /g"; done > brlen.txt    
