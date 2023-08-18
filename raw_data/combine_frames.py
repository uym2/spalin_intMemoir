#! /usr/bin/env python

from sys import argv

inFolder = argv[1]

nframes = 216

for frame in range(1,nframes+1):
    with open(inFolder + "/frame" + str(frame) + ".txt",'r') as fin:
        header = fin.readline().strip().split(",")
        cellID = 1
        for line in fin:
            outline = line.strip() if len(line.strip().split(",")) == 5 else line.strip()+"," 
            print(str(frame) + "," + str(cellID) + "," + outline)
            cellID += 1
