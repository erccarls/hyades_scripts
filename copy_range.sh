#!/bin/bash
        for i in `seq $1 $2`;
        do
                scp "/pfs/carlson/galprop/output/mod_c_$i.hdf5" planck.ucsc.edu:/data/GCE_sys/mod_c/ 
        done    
