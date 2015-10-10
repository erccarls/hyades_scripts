#!/bin/bash
#PBS -N x_co_analysis
#PBS -l nodes=1:ppn=2
#PBS -l pmem=8gb
#PBS -l walltime=8:00:00
#PBS -t 41-48%12
#PBS -q hyper


cd /pfs/carlson/GCE_sys/
python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone.py /pfs/carlson/galprop/output mod_k_"$PBS_ARRAYID" /pfs/carlson/galprop/GALDEF
#python /pfs/carlson/GCE_sys/RunAnalysis.py /pfs/carlson/galprop/output mod_j_"$PBS_ARRAYID"_XCO 0


