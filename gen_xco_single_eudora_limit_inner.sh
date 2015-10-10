#!/bin/bash



for i in {0..10}
do 
	cd /pfs/carlson/GCE_sys/
	python /pfs/carlson/GCE_sys/GenDiffuseModel_X_CO_3zone_P8R2.py /pfs/carlson/galprop/output "$1" /pfs/carlson/galprop/GALDEF $i & 
	#python /pfs/carlson/GCE_sys/RunAnalysis_P8.py /pfs/carlson/galprop/output "$1_$i"_XCO_P8_limit_inner 0 & 
done 


