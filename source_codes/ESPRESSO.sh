#!/bin/sh
#$ -S /bin/sh

# Gao Y, Wang F, Wang R, Kutschera E, Xu Y, Xie S, Wang Y, Kadash-Edmondson KE, Lin L, Xing Y. ESPRESSO: Robust discovery and quantification of transcript isoforms from error-prone long-read RNA-seq data. Sci Adv. 2023 Jan 20;9(3):eabq5072. doi: 10.1126/sciadv.abq5072. Epub 2023 Jan 20. PMID: 36662851; PMCID: PMC9858503.

# only GENCODE38 as guidance
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "ESPRESSO_S.pl, chr${chr}"
perl /path/to/espresso/src/ESPRESSO_S.pl \
-L /path/to/espresso/PB29_merged_chr${chr}_samples.tsv  \
-F /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_chrs/chr${chr}.fna \
-A /path/to/reference/GENCODE/GRCh38/release_38/gtf_chrs/gencode.v38.annotation_chr${chr}.gtf \
-O /path/to/espresso/output_merged_chr${chr} \
-T 5 \
-Q 10

echo "ESPRESSO_C.pl, chr${chr}"
perl /path/to/espresso/src/ESPRESSO_C.pl \
-I /path/to/espresso/output_merged_chr${chr} \
-F /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_chrs/chr${chr}.fna \
-X 0 \
-T 5

echo "ESPRESSO_Q.pl, chr${chr}"
perl /path/to/espresso/src/ESPRESSO_Q.pl \
-L /path/to/espresso/output_merged_chr${chr}/PB29_merged_chr${chr}_samples.tsv.updated \
-A /path/to/reference/GENCODE/GRCh38/release_38/gtf_chrs/gencode.v38.annotation_chr${chr}.gtf \
-V /path/to/espresso/output_merged_chr${chr}/PB29_samples_chr${chr}_N2_R0_compatible_isoform.tsv \
-T 5
done


# GENCODE38 and short-read RNAseq as guidance
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "ESPRESSO_S.pl, chr${chr}"
perl /path/to/espresso/src/ESPRESSO_S.pl \
-L /path/to/espresso/PB29_merged_chr${chr}_samples.tsv  \
-F /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_chrs/chr${chr}.fna \
-A /path/to/reference/GENCODE/GRCh38/release_38/gtf_chrs/gencode.v38.annotation_chr${chr}.gtf \
-B /path/to/sr_RNAseq/all-sr_sum_SJ_thresh3_chr${chr}.out.tab.bed \
-O /path/to/espresso/output_merged_SJ_chr${chr} \
-T 5 \
-Q 10

echo "ESPRESSO_C.pl, chr${chr}"
perl /path/to/espresso/src/ESPRESSO_C.pl \
-I /path/to/espresso/output_merged_SJ_chr${chr} \
-F /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_chrs/chr${chr}.fna \
-X 0 \
-T 5

echo "ESPRESSO_Q.pl, chr${chr}"
perl /path/to/espresso/src/ESPRESSO_Q.pl \
-L /path/to/espresso/output_merged_SJ_chr${chr}/PB29_merged_SJ_chr${chr}_samples.tsv.updated \
-A /path/to/reference/GENCODE/GRCh38/release_38/gtf_chrs/gencode.v38.annotation_chr${chr}.gtf \
-V /path/to/espresso/output_merged_SJ_chr${chr}/PB29_samples_chr${chr}_N2_R0_compatible_isoform.tsv \
-T 5
done 
