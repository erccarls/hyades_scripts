#!/bin/bash
for i in `seq 0 12`; 
do 
python RunAnalysis.py /pfs/carlson/galprop/output/ mod_d_"$i" 0 &
done
