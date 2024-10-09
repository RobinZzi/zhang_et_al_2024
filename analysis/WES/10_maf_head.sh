cd #PATH#/gatk_result/
head -1 sample.hg38_multianno.txt|sed 's/Otherinfo/Tumor_Sample_Barcode/' >header
cat header *maf > patient.maf