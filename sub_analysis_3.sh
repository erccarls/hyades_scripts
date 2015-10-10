#!/bin/bash
#PBS -N ll_analy_3
#PBS -l nodes=1:ppn=1
#PBS -l pmem=3000mb
#PBS -l walltime=16:00:00
#PBS -t 0-70%204


cd /home/carlson/pfs/GCE_sys/
echo "Greetings from job $PBS_ARRAYID"

python RunAnalysis.py /pfs/carlson/galprop/output/ mod_n_"$PBS_ARRAYID"_XCO 3


