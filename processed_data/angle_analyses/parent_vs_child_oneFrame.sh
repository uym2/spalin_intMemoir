grep True ../perFrame_analyses/all_diffussion.txt | awk '{print $1,$2,$6,$7,$3,$4,$5;}' > parent_vs_child_oneFrame.txt
