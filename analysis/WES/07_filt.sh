ref=#PATH#/galtk/Homo_sapiens_assembly38.fasta
bed=#PATH#/hg38.exon.bed

cd #PATH#/gatk_result
for sample in {PT1,PT3,PT4}
do
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  FilterMutectCalls -R $ref\
    -V ${sample}_mutect2.vcf \
    -O ${sample}_somatic.vcf 
    echo "end Mutect2 for ${sample}" `date`
done