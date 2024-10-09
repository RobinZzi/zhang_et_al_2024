bin_trim_galore=trim_galore
dir=#DATAPATH#
cat config|while read id
do
arr=(${id})
fq1=${arr[0]}
fq2=${arr[1]} 
echo  $dir  $fq1 $fq2
nohup $bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 & 
done 
