#! /usr/bin/env python

import os
from problin_libs.sequence_lib import read_sequences
from math import *

variance_t = 3.07 # pre-estimated using all-frame data
t = 215 # tree height is 215 time frames

count_zero = 0
count_all = 0

for d in os.listdir():
    if os.path.isdir(d):
        try:
            msa,_ = read_sequences(d+"/characters.txt",filetype="charMtrx")
            for x in msa:
                count_all += len(msa[x])
                count_zero += len([y for y in msa[x] if y==0])
        except:
            continue        

p = 1-count_zero/count_all
d=-log(count_zero/count_all)
variance_d = t*variance_t/d
print("p="+str(p))
print("d="+str(d))
print("sigma=" + str(sqrt(variance_d)))
