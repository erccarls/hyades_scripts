#!/bin/bash


cd /pfs/carlson/GCE_sys/
python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone_P8R2.py /pfs/carlson/galprop/output $1 /pfs/carlson/galprop/GALDEF
python /pfs/carlson/GCE_sys/RunAnalysis_P8.py /pfs/carlson/galprop/output "$1"_XCO_P8 0
cd /pfs/carlson/galprop/output/
python strip_single.py "$1"_XCO_P8 0 
scp "$1"_XCO_P8_stripped.hdf5 planck:/data/GCE_sys/new


#python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone_P8R2_PSF3.py /pfs/carlson/galprop/output $1 /pfs/carlson/galprop/GALDEF
#python /pfs/carlson/GCE_sys/RunAnalysis_P8_PSF3.py /pfs/carlson/galprop/output "$1"_XCO_P8_PSF3 0 

#python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone.py /pfs/carlson/galprop/output $1 /pfs/carlson/galprop/GALDEF None 1
#python /pfs/carlson/GCE_sys/RunAnalysis.py /pfs/carlson/galprop/output "$1"_XCO 0

#python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone_P8R2_fermi_diffuse.py /pfs/carlson/galprop/output $1 $2




