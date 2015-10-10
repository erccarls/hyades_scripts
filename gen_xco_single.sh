#!/bin/bash


cat <<EOS | qsub -V -q hyper -S /bin/bash -N x_co_analysis_"$1" -l nodes=1:ppn=1,mem=6000mb,walltime=6:00:00 -
cd /pfs/carlson/GCE_sys/
python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone.py /pfs/carlson/galprop/output mod_g_"$1" /pfs/carlson/galprop/GALDEF
python /pfs/carlson/GCE_sys/RunAnalysis.py /pfs/carlson/galprop/output mod_g_"$1"_XCO 0
EOS


