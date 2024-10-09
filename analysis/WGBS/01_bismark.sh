
sample_num=
folder=
#bismark_bam

cd #PATH#/$folder/fastq
qcdir=#PATH#/$folder/qc
tmp_dir=$(mktemp -d)

path=$1
echo $path
files=$(ls $path)
for filename in $files
do
 cd $filename
 (fastqc -t 12 -o $qcdir *.fastq.gz 
 touch "$tmp_dir/$filename.done") &
 cd ..
done

echo "#####################################"
echo "Current date and time: $(date)"
echo "All qc processes are in processing."
echo "#####################################"
       
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All qc processes are completed."
echo "#####################################"

cd #PATH#/$folder/fastq
bin_trim_galore=trim_galore
dir=#PATH#/$folder/trim_result/
tmp_dir=$(mktemp -d)
path=$1
files=$(ls$path)

for filename in $files
do
 cd $filename
 ls *fastq.gz|while read id
 do
 fq1=${filename}_R1.fastq.gz
 fq2=${filename}_R2.fastq.gz
 echo  $fq1 $fq2
 ($bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
  touch "$tmp_dir/$filename.done") &
 done
 cd ..
done

echo "#####################################"
echo "Current date and time: $(date)"
echo "All trim processes are in processing."
echo "#####################################"
       
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All trim processes are completed."
echo "#####################################"

cd #PATH#/bismark/trim_result/${folder}
ref=#PATH#/genome/hg38
dir=#PATH#/bismark/bismark_bam/${folder}
tmp_dir=$(mktemp -d)

i=1
ls *fq.gz |while read id
    do
    headname=${id%_R*}
    if [ $i -eq 2 ];then
    echo ${headname}
    echo $i
    i=1
    (bismark --genome $ref -1 ${headname}_R1.fq.gz -2 ${headname}_R2.fq.gz --o $dir --pbat
    touch "$tmp_dir/$headname.done") &
    else
    i=$[$i+1]
    fi
    done

echo "#####################################"
echo "Current date and time: $(date)"
echo "All bismark_bam are in processing."
echo "#####################################"
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All bismark_bam processes are completed."
echo "#####################################"

#remove_dup
dir=#PATH#/bismark/bismark_bam/${folder}
dir2=#PATH#/bismark/remove_dup/${folder}

cd $dir
tmp_dir=$(mktemp -d)
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
(deduplicate_bismark --bam ${headname}.bam --output_dir $dir2
 touch "$tmp_dir/$headname.done") &
echo $headname
done


echo "#####################################"
echo "Current date and time: $(date)"
echo "All remove_dup are in processing."
echo "#####################################"
while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done
rm -r $tmp_dir

wait
echo "#####################################"
echo "Current date and time: $(date)"
echo "All remove_dup processes are completed."
echo "#####################################"



#meth_ex

dir=#PATH#//$folder/remove_dup
outP=#PATH#//$folder/meth_ex
index=/storage/zhangyanxiaoLab/zhangliwen/genome/hg38


cd $dir
tmp_dir=$(mktemp -d)
ls *.bam |while read id
do
headname=${id%.bam}
echo $headname
(bismark_methylation_extractor -p --no_overlap --parallel 30 --bedGraph --buffer_size 90G --cytosine_report --output ${outP} --genome_folder ${index} ${headname}.bam 
 touch "$tmp_dir/$headname.done") &
done


echo "#####################################"
echo "Current date and time: $(date)"
echo "All meth_ex are in processing."
echo "#####################################"

while [[ $(ls $tmp_dir | wc -l) -lt $sample_num ]];  do
    sleep 5
done

rm -r $tmp_dir

echo "#####################################"
echo "Current date and time: $(date)"
echo "All meth_ex processes are completed."
echo "#####################################"


