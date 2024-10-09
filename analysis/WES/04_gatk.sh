ref=#PATH#/galtk/Homo_sapiens_assembly38.fasta
snp=#PATH#/galtk/dbsnp_146.hg38.vcf
indel=#PATH#/galtk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz



for sample in {NT,PT1,PT2,PT4}

do
	echo $sample

	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
		-I $sample.bam \
		-O ${sample}_marked.bam \
		-M $sample.metrics \
		1>${sample}_log.mark 2>&1

	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
		-I ${sample}_marked.bam \
		-O ${sample}_marked_fixed.bam \
		-SO coordinate \
		1>${sample}_log.fix 2>&1

	#index
	samtools index ${sample}_marked_fixed.bam

	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
		-R $ref  \
		-I ${sample}_marked_fixed.bam  \
		--known-sites $snp \
		--known-sites $indel \
		-O ${sample}_recal.table \
		1>${sample}_log.recal 2>&1

	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./"   ApplyBQSR \
		-R $ref  \
		-I ${sample}_marked_fixed.bam  \
		-bqsr ${sample}_recal.table \
		-O ${sample}_bqsr.bam \
		1>${sample}_log.ApplyBQSR  2>&1

	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
		-R $ref  \
		-I ${sample}_bqsr.bam \
		--dbsnp $snp \
		-O ${sample}_raw.vcf \
		1>${sample}_log.HC 2>&1

done
