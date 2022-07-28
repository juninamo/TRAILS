#!/bin/bash
#$ -S /bin/sh

# download database and script for liftOver
## wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
## mv liftOverPlink.py /path/to/tools/liftOver/

# pre-imputation QC
## remove sample with call rate of <0.98 and all second-degree relations [0.0884], and variants with call rate<0.99, MAF<0.01, HWE P-value≤1.0×10e−5
plink2 --bfile path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685 --king-cutoff 0.0884 --geno 0.01 --maf 0.01 --hwe 1e-5 --set-all-var-ids @_#_\$r_\$a --new-id-max-allele-len 10000 --rm-dup --make-bed --threads 6 --out tmp
plink2 --bfile path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685 --king-cutoff 0.0884 --geno 0.01 --maf 0.01 --hwe 1e-5 --set-all-var-ids @_#_\$r_\$a --new-id-max-allele-len 10000 --rm-dup --make-bed --exclude tmp.rmdup.mismatch --threads 6 --out tmp
plink2 --bfile tmp --mind 0.02 --make-bed --chr 1-22 --threads 6 --out /path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685_QC

for seq_lib in {1..22}; do
plink --bfile /path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685_QC --chr $chr --make-bed --out /path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}
done

# Imputation
## Genetic maps for the 1000 Genomes Project variants; available in https://github.com/joepickrell/1000-genomes-genetic-maps
## 1000GP Phase 3 (n=2,504)
## SHAPEIT → MiniMac3
## Minimac3 Cookbook; https://genome.sph.umich.edu/wiki/Minimac3_Imputation_Cookbook
### pre-phasing
for seq_lib in {1..22}; do
shapeit --input-bed /path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}.bed /path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}.bim /path/to/EGA/plink_files/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}.fam \
--input-map reference/KGP/1000GP_Phase3/genetic_map_chr${seq_lib}_combined_b37.txt \
--output-max /path/to/EGA/Phased/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}_phased.haps /path/to/EGA/Phased/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}_phased.sample \
--thread 8

shapeit -convert \
--input-haps /path/to/EGA/Phased/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}_phased \
--output-vcf /path/to/EGA/Phased/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}_phased_haps.vcf

Minimac3 --refHaps reference/KGP/VCF/G1K_P3/ALL.chr${seq_lib}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz \
--haps /path/to/EGA/Phased/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}_phased_haps.vcf \
--rounds 5 --states 200 \
--format DS,GT,GP \
--prefix /path/to/EGA/Imputed/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}
done

# QC
for i in {1..5}; do
for pop in EUR AFR; do
sed -e "s/^0_//g" /path/to/EGA//${i}_${pop}_sample.list > ${i}_${pop}_sample.list
for seq_lib in {1..22}; do
plink2 \
--vcf /path/to/EGA//GenotypeFiles/Imputed/EvoImmunoPop_Omni5exome_200x3773685_QC_chr${seq_lib}.dose.vcf.gz \
--keep ${i}_${pop}_sample.list \
--set-all-var-ids @_# \
--rm-dup force-first \
--maf 0.05 \
--hwe 1e-6 \
--geno 0.01 \
--minimac3-r2-filter 0.7 \
--make-bed \
--threads 6 \
--allow-extra-chr \
--out /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}
done
done
done

for i in {1..5}; do
for pop in EUR AFR; do
touch mergelist.txt
for seq_lib in {1..22}; do
echo -e "/path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}" >> mergelist.txt
done
plink --merge-list mergelist.txt --make-bed --recode vcf --threads 6 --out /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}
done
done

# liftOver to GRCh38
for i in {1..5}; do
for pop in EUR AFR; do
for seq_lib in {1..22}; do
plink \
--bfile /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib} \
--recode tab \
--threads 6 \
--out /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.tab
/path/to/tools/liftOver/liftOverPlink.py \
-m /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.tab.map \
-p /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.tab.ped \
-o /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38 \
-c /path/to/tools/liftOver/chain_files/hg19ToHg38.over.chain \
-e /path/to/tools/liftOver/liftOver
plink --file /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38 \
--make-bed \
--allow-extra-chr \
--threads 6 \
--out /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38

# correct variant id to chr_pos
plink2 \
--bfile /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38 \
--set-all-var-ids @_# \
--rm-dup force-first \
--maf 0.05 \
--hwe 1e-6 \
--geno 0.01 \
--make-bed \
--threads 6 \
--allow-extra-chr \
--out tmp/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38
mv tmp/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38* /path/to/EGA//GenotypeFiles/Imputed/
done
done
done

# Check how much the variants are
for i in {1..5}; do
for pop in EUR AFR; do
for seq_lib in {1..22}; do
wc -l /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.hg38.bim
wc -l /path/to/EGA//GenotypeFiles/Imputed/QC_P3_R20.7_MAF0.05_HWE1E06_nodup_Mono-${i}_${pop}_chr${seq_lib}.bim
done
done
done
