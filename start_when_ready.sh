

while true; 
do
    num=`qstat -a | grep galprop_mod | wc -l`
    test $num -eq 0 && break
    sleep 1
done

cd /pfs/carlson/GCE_sys/
./sub_analysis.sh
