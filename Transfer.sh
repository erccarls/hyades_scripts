# cd ..//galprop/output/
# scp mod_s_42_XCO_P8.hdf5 planck:/data/GCE_sys/new & 
# scp mod_s_46_XCO_P8.hdf5 planck:/data/GCE_sys/new 
# cd ../../GCE_sys/


cd ..//galprop/output/

# for i in {42..48}
# do
#  	python strip_single.py mod_s_"$i"_XCO_P8_corrected 0
# done

#  for i in {42..48}
#  do
#  	scp mod_s_"$i"_XCO_P8_corrected_stripped.hdf5 planck:/data/GCE_sys/new 
#  done




for i in {0..101}
do
 	python strip_single.py mod_p_"$i"_XCO_P8_corrected 0
done

for i in {0..101}
do
	scp mod_p_"$i"_XCO_P8_corrected_stripped.hdf5 planck:/data/GCE_sys/new 
done






# for i in {28..41}
# do
#  	python strip_single.py mod_s_"$i"_XCO_P8_corrected 0
# done

#  for i in {28..41}
#  do
#  	scp mod_s_"$i"_XCO_P8_corrected_stripped.hdf5 planck:/data/GCE_sys/new 
#  done




# for i in {0..12}
# do
#  	python strip_single.py mod_q_"$i"_XCO_P8_corrected 0
# done

#  for i in {0..12}
#  do
#  	scp mod_q_"$i"_XCO_P8_corrected_stripped.hdf5 planck:/data/GCE_sys/new 
#  done




# for i in {0..10}
# do
# 	python strip_single.py mod_s6_"$i"_XCO_P8_corrected 0 
# 	#python strip_single.py mod_s5_"$i"_XCO_P8_corrected_freeH 0 
# done

# for i in {0..10}
# do
# 	scp mod_s6_"$i"_XCO_P8_corrected_stripped.hdf5 planck:/data/GCE_sys/new 
# 	#scp mod_s6_"$i"_XCO_P8_corrected_freeH_stripped.hdf5 planck:/data/GCE_sys/new 
# done



#scp mod_s_44_XCO_P8_corrected.hdf5 planck:/data/GCE_sys/new & 
#scp mod_s_47_XCO_P8_corrected.hdf5 planck:/data/GCE_sys/new 
cd ../../GCE_sys/

#cd ..//galprop/output/
# scp mod_sMS04_42_XCO_P8.hdf5 planck:/data/GCE_sys/new & 
# scp mod_sMS04_46_XCO_P8.hdf5 planck:/data/GCE_sys/new 

#for i_model in {43,47,48}
#do
#	scp mod_sMS04_"$i_model"_XCO_P8.hdf5 planck:/data/GCE_sys/new 
#	scp mod_sMS04_"$i_model"_XCO_P8.hdf5 planck:/data/GCE_sys/new 
#done 
#cd ../../GCE_sys/

