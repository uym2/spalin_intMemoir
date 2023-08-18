#!/bin/bash

for x in *_frames ; do
    y=`echo $x | sed -e "s/_frames//g"`
    python combine_frames.py $x | sed -e "s/^/$y,/g"
done    
