#! /bin/bash

smooth=$1

python infer_rates.py problin_k10_tree_const.nwk $smooth | tail -n+3 | sed -e "s/^/constant /g" >> problin_k10_inferred_rates_smooth_$smooth.txt

python infer_rates.py problin_k10_tree_lin.nwk $smooth | tail -n+3 | sed -e "s/^/linear /g" >> problin_k10_inferred_rates_smooth_$smooth.txt

python infer_rates.py problin_k10_tree_exp.nwk $smooth | tail -n+3 | sed -e "s/^/exp /g" >> problin_k10_inferred_rates_smooth_$smooth.txt

python infer_rates.py problin_k10_tree_exp2.nwk $smooth | tail -n+3 | sed -e "s/^/exp2 /g" >> problin_k10_inferred_rates_smooth_$smooth.txt
