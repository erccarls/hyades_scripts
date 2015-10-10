#!/bin/bash
#PBS -N ll_analy_10
#PBS -l nodes=1:ppn=2
#PBS -l pmem=6000mb
#PBS -q hyper
#PBS -l walltime=12:00:00
#PBS -t 0-56%57


cd /home/carlson/pfs/GCE_sys/
echo "Greetings from job $PBS_ARRAYID"

python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$PBS_ARRAYID"_XCO_P8 0

