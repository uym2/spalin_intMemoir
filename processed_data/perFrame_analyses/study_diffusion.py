#! /usr/bin/env python
from math import sqrt

infile = "all_tracked.txt"

nframes = 216
D = {} # will be a dictionary of lists of dictionaries

with open(infile,'r') as fin:
    fin.readline() # skip the header
    for line in fin:
        ID,frame,cell,pos_x,pos_y,parent,child_1,child_2 = line.strip().split(",")
        if ID not in D:
            D[ID] = []
            for f in range(nframes):
                D[ID].append({})
        #D[ID][int(frame)-1]  = {}   
        D[ID][int(frame)-1][cell] = (float(pos_x),float(pos_y),parent,child_1,child_2)

for ID in D:
    for frame in range(1,nframes):
        for cell in D[ID][frame]:
            cell_par = D[ID][frame][cell][2] # get the parent cell
            is_div =  D[ID][frame-1][cell_par][4] != ''
            print(ID,frame+1,cell,D[ID][frame][cell][0],D[ID][frame][cell][1],D[ID][frame-1][cell_par][0],D[ID][frame-1][cell_par][1],is_div)
