#! /bin/bash

smooth=$1

python infer_rates.py true_scaled_tree_const.nwk $smooth | tail -n+3 | sed -e "s/^/constant /g" > true_scaled_inferred_rates_smooth_$smooth.txt

python infer_rates.py true_scaled_tree_lin.nwk $smooth | tail -n+3 | sed -e "s/^/linear /g" >> true_scaled_inferred_rates_smooth_$smooth.txt

python infer_rates.py true_scaled_tree_exp.nwk $smooth | tail -n+3 | sed -e "s/^/exp /g" >> true_scaled_inferred_rates_smooth_$smooth.txt

python infer_rates.py true_scaled_tree_exp2.nwk $smooth | tail -n+3 | sed -e "s/^/exp2 /g" >> true_scaled_inferred_rates_smooth_$smooth.txt
