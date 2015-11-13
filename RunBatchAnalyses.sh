# for i_model in {42,46}
# do
# 	#for i_analysis in {0,15,16,19,20,21}
# 	for i_analysis in {19,}
# 	do
# 		for i_psf in {-1..3}
# 		do	
# 			while true; 
# 				do
# 				   num=`ps -u carlson | grep python | wc -l`
# 				   test $num -le 24 && break
# 				   sleep 1
# 				done
# 			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8 "$i_analysis" PSF="$i_psf" & 

# 		done	
# 	done
# done




# for i_model in {42,46}
# do
# 	#for i_analysis in {0,15,16,19,20,21}
# 	for i_analysis in {0,}
# 	do
# 		for i_psf in {-1..3}
# 		do	
# 			while true; 
# 				do
# 				   num=`ps -u carlson | grep python | wc -l`
# 				   test $num -le 24 && break
# 				   sleep 1
# 				done
# 			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8 "$i_analysis" PSF="$i_psf" & 

# 		done	
# 	done
# done


# Transfer Files.
# while true; 
# do
#    num=`ps -u carlson | grep RunAnalysis | wc -l`
#    test $num -le 0 && break
#    sleep 1
# done
# cd ../galprop/output/
# scp mod_s_42_XCO_P8.hdf5 planck:/data/GCE_sys/new & 
# scp mod_s_46_XCO_P8.hdf5 planck:/data/GCE_sys/new
# cd ../../GCE_sys/





# ------------------------#------------------------#------------------------
# Hyades runs

#for i_model in {28..55}
# for i_model in {42..48}
# do
# 	for i_analysis in {0,}
# 	#for i_analysis in {3,5,7,8,10,}
# 	#for i_analysis in {22,23}
# 	do
# 		for i_psf in {-1,}
# 		do	
# 			cat <<EOS | qsub -V -q normal -S /bin/bash -N analysis_"$i_analysis" -l nodes=1:ppn=1,walltime=24:00:00 - 
# 			export OMP_NUM_THREADS=1
# 			cd /pfs/carlson/GCE_sys/
# 			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8_corrected "$i_analysis" 0 PSF="$i_psf" 
# 			#python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8_corrected 0 0 PSF=-3 
# EOS

# 		done
# 	done
# done



#---------------------------
# Hyades mod_n!!!!

for i_model in {42..48}
do
	for i_analysis in {0, 24,}
	#for i_analysis in {3,5,7,8,10,}
	#for i_analysis in {22,23}
	do
		for i_psf in {-1,}
		do	
			cat <<EOS | qsub -V -q normal -S /bin/bash -N analysis_"$i_analysis" -l nodes=1:ppn=1,walltime=6:00:00 - 
			export OMP_NUM_THREADS=1
			cd /pfs/carlson/GCE_sys/
			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8_corrected_freeH "$i_analysis" 0 PSF="$i_psf" 
EOS

		done
	done
done




# while true; 
# do
#    num=`qstat -u carlson | grep analysis | wc -l`
#    test $num -le 0 && break
#    sleep 1
# done

# echo "Done Running Analysis" | mail -s "Done Running Analysis" 5034421563@vtext.com
# #------------------------#------------------------#------------------------#------------------------







#------------------------#------------------------#------------------------
# Hyades runs

#for i_model in {14,15,16,17,18,19,20,42,43,44,45,46,47,48}

# for i_model in {42,43,44,45,46,47,48}
# do
# 	for i_analysis in {22,23}
# 	#for i_analysis in {1,2,3,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21}
# 	do
# 		for i_psf in {-1,3}
# 		do	

# 			while true; 
# 			do
# 			   num=`qstat -u carlson | grep analysis_ | wc -l`
# 			   test $num -le 128 && break
# 			   sleep 1
# 			done

# 			cat <<EOS | qsub -V -q normal -S /bin/bash -N analysis_"$i_analysis" -l nodes=1:ppn=1,walltime=16:00:00 - 
# 			export OMP_NUM_THREADS=1
# 			cd /pfs/carlson/GCE_sys/
# 			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_s_"$i_model"_XCO_P8_corrected "$i_analysis" PSF="$i_psf" 
# EOS

# 		done
# 	done
# done

# while true; 
# do
#    num=`qstat -u carlson | grep analysis | wc -l`
#    test $num -le 0 && break
#    sleep 1
# done

# echo "Done Running Analysis for corrected models" | mail -s "Done Running Analysis" erccarls@ucsc.edu





# for i_model in {43,47,48}
# do
# 	for i_analysis in {0,15,16,19,20,21}
# 	#for i_analysis in {21,}
# 	do
# 		for i_psf in {-1..3}
# 		do	
# 			cat <<EOS | qsub -V -q normal -S /bin/bash -N analysis_"$i_analysis" -l nodes=1:ppn=1,walltime=1:00:00 - 
# 			export OMP_NUM_THREADS=1
# 			cd /pfs/carlson/GCE_sys/
# 			python RunAnalysis_P8.py /pfs/carlson/galprop/output/ mod_sMS04_"$i_model"_XCO_P8 "$i_analysis" PSF="$i_psf" 
# EOS

# 		done
# 	done
# done

# while true; 
# do
#    num=`qstat -u carlson | grep analysis | wc -l`
#    test $num -le 0 && break
#    sleep 1
# done

#echo "Done Running Analysis for corrected models" | mail -s "Done Running Analysis" 5034421563@vtext.com





