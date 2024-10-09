GENOME=#PATH#/genome/galtk/Homo_sapiens_assembly38.fasta
bed=#PATH#/genome/galtk/hg38.exon.bed
cd #PATH#/gatk_result/
cnvkit.py batch PT*bqsr.bam \
--normal  NT_bqsr.bam \
--targets  ${bed} \
--fasta $GENOME  \
--drop-low-coverage --scatter --diagram --method amplicon \
--output-reference my_reference.cnn --output-dir #PATH#/cnvkit_result
