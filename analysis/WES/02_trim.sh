bin_trim_galore=trim_galore
dir=#DATAPATH#/trim_result/
i=1
ls *fq.gz|while read id
do
if [ $i -eq 2 ];then
    i=1
else
headname=${id%_R*}
fq1=${headname}_R1.fq.gz
fq2=${headname}_R2.fq.gz
echo  $dir  $fq1 $fq2
$bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 & 
i=$[$i+1]
fi
done 

