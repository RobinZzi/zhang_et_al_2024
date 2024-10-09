ref=#PATH#/galtk/Homo_sapiens_assembly38.fasta
bed=#PATH#/genome/galtk/hg38.exon.bed

cd #PATH#/gatk_result
    gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R $ref \
    -I sample_bqsr.bam -tumor sample \
    -I NT_bqsr.bam -normal NT \
    -L $bed  \
    -O sample_mutect2.vcf


    
    
    
