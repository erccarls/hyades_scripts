#!/bin/bash


cat <<EOS | qsub -q normal -S /bin/bash -N xco_"$1" -l nodes=1:ppn=1,mem=10gb,walltime=6:00:00 -
cd /pfs/carlson/GCE_sys/
python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone_P8R2.py /pfs/carlson/galprop/output $1 /pfs/carlson/galprop/GALDEF
#python /pfs/carlson/GCE_sys/RunAnalysis_P8.py /pfs/carlson/galprop/output "$1"_XCO_P8 0
EOS


