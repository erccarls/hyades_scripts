#!/bin/bash
#PBS -N ll_analy_1
#PBS -l nodes=1:ppn=1
#PBS -l pmem=3000mb
#PBS -l walltime=12:00:00
#PBS -t 0-70%204


cd /home/carlson/pfs/GCE_sys/
echo "Greetings from job $PBS_ARRAYID"

python RunAnalysis.py /pfs/carlson/galprop/output/ mod_m_"$PBS_ARRAYID"_XCO 1


