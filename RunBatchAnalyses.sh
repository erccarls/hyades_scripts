for i_model in {42,46}
do
	#for i_analysis in {0,15,16,19,20,21}
	for i_analysis in {19,}
	do
		for i_psf in {-1..3}
		do	
			while true; 
				do
				   num=`ps -u carlson | grep python | wc -l`
				   test $num -le 24 && break
				   sleep 1
				done
			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8 "$i_analysis" PSF="$i_psf" & 

		done	
	done
done