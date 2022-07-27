#!/bin/sh
#$ -S /bin/sh

# flair

PBs=(
    "PB01"
    "PB02"
    "PB03"
    "PB04"
    "PB05"
    "PB06"
    "PB07"
    "PB08"
    "PB09"
    "PB10"
    "PB11"
    "PB12"
    "PB13"
    "PB14"
    "PB15"
    "PB16"
    "PB17"
    "PB18"
    "PB19"
    "PB20"
    "PB21_A"
    "PB21_B"
    "PB21_C"
    "PB22"
    "PB23"
    "PB24"
    "PB25"
    "PB26"
    "PB27"
    "PB28"
    "PB29"
)

cells=(
    "NaiveCD4"
	"Th1"
	"Th2"
	"Th17"
	"Tfh"
	"Fra1"
	"Fra2-aTreg"
	"Fra3"
	"LAG3Treg"
	"MemoryCD4"
	"Thx"
	"NaiveCD8"
	"CD8effector"
	"CD8centralmem"
	"CD8effectormem"
	"NaiveB"
	"unswmemoryB"
	"swmemoryB"
	"DNB"
	"plasmablast"
	"plasmacytoidDC_A"
	"plasmacytoidDC_B"
	"plasmacytoidDC_C"
	"myeloidDC"
	"NK"
	"monocyteCD16"
	"monocyteCD16minus"
	"nonclassicalMonocyte"
	"intermediateMonocyte"
	"PBMC"
	"Neutrophil"
)

for n in {0..30} ; do
echo "merge, ${PBs[n]}; ${cells[n]}"
mkdir --parents /path/to/${cells[n]}/fasta/
zcat /path/to/${PBs[n]}*/pass/fastq_*.fastq.gz > /path/to/${cells[n]}/fasta/${cells[n]}.fq
done


for n in {0..30} ; do
mkdir --parents /path/to/${cells[n]}/flair
done

# alighn
for n in {0..30} ; do
echo "alighn GRCh38, ${PBs[n]}; ${cells[n]}"
cd /path/to/${cells[n]}/flair/
/path/to/bio/python3/bin/python3 /path/to/tools/flair/flair.py align \
--genome /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi \
--reads /path/to/${cells[n]}/fasta/${cells[n]}.fq \
--minimap2 /path/to/tools/minimap2/minimap2 \
--samtools /opt/bin/samtools \
--quality 10 \
--threads 30 \
--output GRCh38.aligned
done

# correct alighned.bed by junction reads from short-read RNA-seq
for n in {0..30} ; do
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "correct, ${PBs[n]}; ${cells[n]}, chr${chr}"
cd /path/to/${cells[n]}/flair/
/path/to/bio/python3/bin/python3 /path/to/tools/flair/flair.py correct \
--genome /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
--query GRCh38.aligned.bed \
--shortread /path/to/sr_RNAseq/all-sr_sum_SJ_thresh3_chr${chr}.out.tab \
--ss_window 10 \
--threads 30 \
--print_check \
--output GRCh38_all-sr-correct_chr${chr}
done
done > /path/to/`date "+%Y%m%d"`_flair_all-sr-correct.log 2>&1


# merge all_corrected.bed from all cell-types
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "chr${chr}"
LANG=C sort -V -k 1,1 -k 2,2n /path/to/*/flair/GRCh38_all-sr-correct_chr${chr}_all_corrected.bed > /path/to/PB29_common/chrs/GRCh38_all-sr-correct_chr${chr}_all_corrected.bed
done


# merge refTSS and TSS_classifier(>relaxed) and split by chrosome
cat \
/path/to/reference/refTSS/v3.3/human/refTSS_v3.3_human_coordinate.hg38.bed \
/path/to/reference/FANTOM5/TSS_classifier/TSS_human_hg38_relaxed_strict_TSS.bed \
| LANG=C sort -V -k 1,1 -k 2,2n \
| bedtools merge -i - > /path/to/annotation/TSS/merged_refTSS-v3.3_TSS-classifier-relaxed-strict_hg38.bed

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "chr${chr}"
fgrep -w "chr${chr}" /path/to/annotation/TSS/merged_refTSS-v3.3_TSS-classifier-relaxed-strict_hg38.bed > /path/to/annotation/TSS/chrs/merged_refTSS-v3.3_TSS-classifier-relaxed-strict_hg38_chr${chr}.bed
done


# collapse corrected.bed
cd /path/to/PB29_common/
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "collapse chr${chr}"
/path/to/bio/python3/bin/python3 /path/to/tools/flair/flair.py collapse \
--genome /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set_chrs/chr${chr}.fna \
--reads /path/to/*/fasta/*.fq \
--query /path/to/PB29_common/chrs/GRCh38_all-sr-correct_chr${chr}_all_corrected.bed \
--gtf /path/to/reference/GENCODE/GRCh38/release_38/gtf_chrs/gencode.v38.annotation_chr${chr}.gtf \
--minimap2 /path/to/tools/minimap2/minimap2 \
--mm2_args=-I8g,--MD \
--samtools /opt/bin/samtools \
--bedtools /opt/bin/bedtools \
--promoters /path/to/annotation/TSS/chrs/merged_refTSS-v3.3_TSS-classifier-relaxed-strict_hg38_chr${chr}.bed \
--end_window 100 \
--support 3 \
--stringent \
--no_redundant none \
--max_ends 2 \
--filter default \
--threads 30 \
--quality 10 \
--keep_intermediate \
--temp_dir ./chr${chr} \
--output GRCh38_all-sr-correct_chr${chr}
done > /path/to/`date "+%Y%m%d"`_flair_all-sr-correct_collapse_common.log 2>&1

# merge all chrosomes
cat /path/to/PB29_common/GRCh38_all-sr-correct_chr*.isoforms.gtf | awk '$4 < $5{print $0}' > /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.gtf
cat /path/to/PB29_common/GRCh38_all-sr-correct_chr*.isoforms.fa > /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.fa
/path/to/tools/minimap2/minimap2 -d /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.mmi /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.fa
cat /path/to/PB29_common/GRCh38_all-sr-correct_chr*.isoforms.bed | awk '$2 < $3{print $0}' | grep -v ",-" > /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.bed


# quantify

touch /path/to/reads_manifest.tsv
for n in {0..30} ; do
paste \
<(echo -e ${cells[n]}) \
<(echo -e $(grep -w ${cells[n]} /path/to/reads_manifest_wildcard.tsv | cut -f 2)) \
| paste - <(echo -e "batch1") \
| paste - <(echo -e "/path/to/${cells[n]}/fasta/${cells[n]}.fq") >> /path/to/reads_manifest.tsv
done

cd /path/to/PB29_common/

/path/to/bio/python3/bin/python3 /path/to/tools/flair/flair.py quantify \
--reads_manifest /path/to/reads_manifest.tsv \
--isoforms /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.fa \
--minimap2 /path/to/tools/minimap2/minimap2 \
--samtools /opt/bin/samtools \
--quality 10 \
--tpm \
--threads 30 \
--output GRCh38_all-sr-correct_counts_matrix.tsv \
> /path/to/`date "+%Y%m%d"`_flair_quantify_GRCh38_all-sr-correct_common.log 2>&1

# predictProductivity
/path/to/bio/python3/bin/python3 /path/to/tools/flair/bin/predictProductivity_modified.py \
--input_isoforms /path/to/PB29_common/GRCh38_all-sr-correct.isoforms.bed \
--gtf /path/to/reference/GENCODE/GRCh38/release_38/gencode.v38.annotation.gtf \
--genome_fasta /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
--longestORF --append_column > /path/to/PB29_common/GRCh38_all-sr-correct_productivity.bed

# mark_intron_retention
/path/to/bio/python3/bin/python3 /path/to/tools/flair/bin/mark_intron_retention.py \
/path/to/PB29_common/GRCh38_all-sr-correct.isoforms.bed \
/path/to/PB29_common/GRCh38_all-sr-correct.isoforms.mark_intron_retention.psl \
/path/to/PB29_common/GRCh38_all-sr-correct.isoforms.mark_intron_retention_coords.txt




