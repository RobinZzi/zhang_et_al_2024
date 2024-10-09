ref=#PATH#/genome/hg38.fa
output_dir=#PATH#/bwa_result
i=1
ls *gz  |while read id
do
    if [ $i -eq 2 ];then
    i=1
    else
    headname=${id%_R*}
    fq1=${headname}_R1_val_1.fq.gz
    fq2=${headname}_R2_val_2.fq.gz
    echo $fq1
    echo $fq2
    bwa mem -M -t 8 -R "@RG\tID:$headname\tSM:$headname\tLB:WES\tPL:Illumina" $ref $fq1 $fq2 | samtools sort -@ 11 -m 1G --reference $fa  -o ${output_dir}/${headname}.bam &
    i=$[$i+1]
    fi
done




