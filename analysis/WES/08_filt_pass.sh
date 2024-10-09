cd #PATH#/gatk_result
for i in *_somatic.vcf
do
j=$(basename "$i" _somatic.vcf )
echo $j
cat $i | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ${j}_filter.vcf
done