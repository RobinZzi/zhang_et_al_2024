annovar=#PATH#/annovar/table_annovar.pl
ref=#PATH#/genome/humandb/
cd #PATH#/gatk_result/
for id in {PT1,PT2,PT3,PT4}
do
$annovar ${id}_filter.vcf $ref \
-buildver hg38 \
-out ${id} \
-remove \
-protocol refGene,knownGene,clinvar_20170905 \
-operation g,g,f \
-nastring . \
-vcfinput
done