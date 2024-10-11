cd #PATH#/data/hcc/hcc1/10x

#count
sample="HCC1-10X-NT"

cellranger count --id=${sample} \
                    --transcriptome=#PATH#/database/cellrgindex/refdata-gex-GRCh38-2020-A \
                    --fastqs=#PATH#/data/hcc/hcc1/10x \
                    --sample=${sample}-1,${sample}-2,${sample}-3,${sample}-4 \
                    --localcores=30 \
                    --disable-ui


#aggr
sample="HCC1-10X"
cellranger aggr --id=${sample} \
                   --csv=#PATH#/data/hcc/hcc1/10x/aggr.lib.csv \
                   --normalize=mapped \
                   --localcores=30 \
                   --disable-ui