#!/bin/bash
#$ -S /bin/sh

# download from 1000 Genomes on GRCh38 and convert ID of GEUVADIS
## https://www.internationalgenome.org/data-portal/data-collection/grch38

# QC, merge and get AF
for pop in EUR YRI; do
for seq_lib in {1..22}; do
plink2 \
--vcf /path/to/GEUV/VCF/GEUVADIS_GRCh38_${pop}_chr${seq_lib}.vcf.gz \
--set-all-var-ids @_# \
--rm-dup force-first \
--maf 0.05 \
--hwe 1e-6 \
--geno 0.01 \
--make-bed \
--threads 6 \
--out /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}_chr${seq_lib}
plink2 --bfile /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}_chr${seq_lib} \
--freq \
--out /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}_chr${seq_lib}
done

touch mergelist.txt
for seq_lib in {1..22}; do
echo -e "/data03/inamo//path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}_chr${seq_lib}" >> mergelist.txt
done

plink --merge-list mergelist.txt --make-bed --recode vcf --threads 6 --out /data03/inamo//path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop} --allow-extra-chr

bgzip -f /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}.vcf && tabix -p vcf /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}.vcf.gz
/path/to/QTLtools pca --vcf /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}.vcf.gz --out /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop} --center --scale --maf 0.05 --distance 5000
head /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}.pca -n 11 > /path/to/GEUV/VCF/QC_MAF0.05_HWE1E06_nodup_GRCh38_${pop}_PC10.txt

done