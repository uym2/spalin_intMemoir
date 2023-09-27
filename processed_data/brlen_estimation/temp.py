from problin_libs.sequence_lib import read_sequences
from sys import argv
from math import *

data = argv[1] 
def count_missing(mtx):
    nzeros, nmissing, total = 0, 0, 0 
    for c in mtx:
        seq = mtx[c]
        nzeros += sum([1 if (ch == 0 or ch == '0') else 0 for ch in seq])
        nmissing += sum([1 if ch == '?' else 0 for ch in seq])
        total += len(seq)
    return nzeros, nmissing, total

msa,_ = read_sequences(data+"/characters.txt",filetype="charMtrx")
nzeros,nmissing,total = count_missing(msa)
print(data,total,-log(nzeros/total))
