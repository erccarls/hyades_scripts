#!/bin/bash
#PBS -N ll_analy_4
#PBS -l nodes=1:ppn=1
#PBS -l pmem=3000mb
#PBS -l walltime=48:00:00
#PBS -t 0-129%200


cd /home/carlson/pfs/GCE_sys/
echo "Greetings from job $PBS_ARRAYID"

python RunAnalysis.py /pfs/carlson/galprop/output/ mod_k_"$PBS_ARRAYID"_XCO 4


