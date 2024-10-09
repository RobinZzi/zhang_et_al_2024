data_path=#PATH#/bismark/sorted_bam
output_path=#PATH#/bismark/samtools_depth_result
folder=
bed=#PATH#/genome/galtk/hg38.exon.bed

cd ${data_path}/${folder}
for bam in *.bam; do 
samtools depth -b ${bed} ${bam} > ${output_path}/${folder}/${bam%.bam}.depth
done

echo '##############################'
echo "Current date and time: $(date)"
echo '############done##############'
echo '##############################'