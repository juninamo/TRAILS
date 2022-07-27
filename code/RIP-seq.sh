#!/bin/bash
#$ -S /bin/sh

# LCL (SunyAlbany_RipSeq_GM12878_ELAVL1)
## GSE35585; GSM944520

#============================================
# step 0: download
#============================================

GEO_id="GSM944520"
echo $GEO_id

for rna_id in SRR504447 SRR504448; do
mkdir -p /path/to/RipSeq/${GEO_id}/dump/${rna_id}
if [ -e "/path/to/RipSeq/${GEO_id}/dump/${rna_id}/${rna_id}.fastq.gz" ]; then
echo -e "${rna_id}.fastq.gz exists"
else
echo -e "${rna_id}.fastq.gz not exists"
DIR="/path/to/RipSeq/${GEO_id}/dump/${rna_id}"
tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump-orig.2.11.2 -O ${DIR} --split-3 -e 4 -p ${rna_id}
fi
done
cat /path/to/RipSeq/${GEO_id}/dump/*/*fastq > /path/to/RipSeq/${GEO_id}/dump/merged.fastq
rm /path/to/RipSeq/${GEO_id}/dump/*/*fastq
tools/pigz-2.7/pigz -p 20 /path/to/RipSeq/${GEO_id}/dump/merged.fastq

eval "$(/data02/home/juninamo/miniconda3/bin/conda shell.bash hook)"
# conda create -n fastqc
conda activate fastqc
# conda install -c bioconda fastqc

GEO_id="GSM944520"
mkdir -p /path/to/RipSeq/${GEO_id}/fastqc
fastqc -t 8 -o /path/to/RipSeq/${GEO_id}/fastqc /path/to/RipSeq/${GEO_id}/dump/merged.fastq.gz -f fastq

trim_galore -q 30 -o /path/to/RipSeq/${GEO_id}/trim_galore/ /path/to/RipSeq/${GEO_id}/dump/merged.fastq.gz

# remove rRNA + Repetitive elements
## https://imamachi-n.hatenablog.com/entry/2017/03/31/215218
## cat reference/fasta/contam_Ribosomal_RNA.fa /data02/home/juninamo/reference/RepeatMasker/RepBaseRepeatMaskerEdition-20181026/Libraries/RepBase.fasta > reference/fasta/contam_RepBase_Ribosomal_RNA.fa
## mkdir -p tools/STAR/RepBase_rRNA_contam
## /path/to/STAR \
## --runThreadN 8 \
## --runMode genomeGenerate \
## --genomeDir tools/STAR/RepBase_rRNA_contam \
## --limitGenomeGenerateRAM 40000000000 \
## --genomeFastaFiles reference/fasta/contam_RepBase_Ribosomal_RNA.fa
mkdir -p /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna
/path/to/STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir tools/STAR/RepBase_rRNA_contam \
--readFilesIn /path/to/RipSeq/${GEO_id}/trim_galore/merged_trimmed.fq.gz \
--readFilesCommand gunzip -c \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/rm_repbase_rrna.fastq \
--outSAMattributes All \
--outStd BAM_Unsorted \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--alignEndsType EndToEnd > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/repbase_rrna_comtam.bam
mv /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/rm_repbase_rrna.fastqUnmapped.out.mate1 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/rm_repbase_rrna.fastq
/path/to/STAR \
--runThreadN 10 \
--genomeDir tools/STAR/genome_GRCh38_PB29/ \
--readFilesIn /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/rm_repbase_rrna.fastq \
--outFileNamePrefix /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--outFilterType BySJout \
--outFilterScoreMin 10 \
--sjdbGTFfile /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
--alignEndsType EndToEnd
/path/to/samtools sort -@ 8 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/Aligned.sortedByCoord.out.bam -o /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bam
/path/to/samtools index -@ 8 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bam
rm /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/Aligned.sortedByCoord.out.bam
/path/to/samtools view -h /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.sam
/path/to/samtools view -@ 10 -q 255 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bam -o /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam
/path/to/samtools index -@ 8 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam
/path/to/samtools view -h /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.sam


# QC
picard QualityScoreDistribution INPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.bam OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.quality_score_dist.tsv CHART_OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.quality_score_dist.pdf
picard QualityScoreDistribution INPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_uniq.bam OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_uniq.quality_score_dist.tsv CHART_OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_uniq.quality_score_dist.pdf
picard QualityScoreDistribution INPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.bam OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.quality_score_dist.tsv CHART_OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.quality_score_dist.pdf
picard QualityScoreDistribution INPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.bam OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.quality_score_dist.tsv CHART_OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.quality_score_dist.pdf
picard QualityScoreDistribution INPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bam OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.quality_score_dist.tsv CHART_OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.quality_score_dist.pdf
picard QualityScoreDistribution INPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.quality_score_dist.tsv CHART_OUTPUT=/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.quality_score_dist.pdf
bamtools stats -in /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.bamtools.stats
bamtools stats -in /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_uniq.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_uniq.bamtools.stats
bamtools stats -in /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.bamtools.stats
bamtools stats -in /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.bamtools.stats
bamtools stats -in /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_sorted.bamtools.stats
bamtools stats -in /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bamtools.stats

# Number of reads mapped to each position in each transcript
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.depth.tab
# Number of reads mapped to each transcript
/path/to/samtools idxstats -@ 8 /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.bam > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.idxstats.tab

# Visualization (For UCSC genome browser)
bedtools genomecov -ibam /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.bam -bg -split > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.bg
bedtools genomecov -ibam /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.bam -bg -split > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.bg
bedtools genomecov -ibam /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam -bg -split > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bg
echo "track type=bedGraph name=LCL_ELAVL1_PAR-CLIP_sorted description=LCL_ELAVL1_PAR-CLIP_sorted visibility=2 maxHeightPixels=40:40:20 color=128,0,128" > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted_for_UCSC.txt
echo "track type=bedGraph name=LCL_ELAVL1_PAR-CLIP_uniq description=LCL_ELAVL1_PAR-CLIP_uniq visibility=2 maxHeightPixels=40:40:20 color=128,0,128" > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq_for_UCSC.txt
echo "track type=bedGraph name=LCL_ELAVL1_PAR-CLIP_rm_repbase_rrna_uniq description=LCL_ELAVL1_PAR-CLIP_rm_repbase_rrna_uniq visibility=2 maxHeightPixels=40:40:20 color=128,0,128" > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq_for_UCSC.txt
cat /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted_for_UCSC.txt /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted.bg > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted_for_UCSC.bg
cat /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq_for_UCSC.txt /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq.bg > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq_for_UCSC.bg
cat /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq_for_UCSC.txt /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bg > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq_for_UCSC.bg
bzip2 -c /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted_for_UCSC.bg > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_sorted_for_UCSC.bg.bz2
bzip2 -c /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq_for_UCSC.bg > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome/genome_uniq_for_UCSC.bg.bz2
bzip2 -c /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq_for_UCSC.bg > /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq_for_UCSC.bg.bz2

# rf-count
tools/RNAFramework/rf-count \
-p 30 -wt 30 \
-s /path/to/samtools -r \
-o /path/to/RipSeq/${GEO_id}/RNA_Framework/count_GEUVexp_nonredundunt -ow \
-f /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.filtered.fasta \
/path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29/transcriptome_sorted.bam

# rf-rctools
tools/RNAFramework/rf-rctools view /path/to/RipSeq/${GEO_id}/RNA_Framework/count_GEUVexp_nonredundunt/transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/RipSeq/${GEO_id}/RNA_Framework/count_GEUVexp_nonredundunt/transcriptome_sorted.tab
tools/RNAFramework/rf-rctools stats /path/to/RipSeq/${GEO_id}/RNA_Framework/count_GEUVexp_nonredundunt/transcriptome_sorted.rc > /path/to/RipSeq/${GEO_id}/RNA_Framework/count_GEUVexp_nonredundunt/transcriptome_sorted.stats

### Extract 3UTR bed file
## https://imamachi-n.hatenablog.com/entry/2017/03/29/211905
perl tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/bed12to3UTRbed.pl \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.bed \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_3UTR.bed \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_non-3UTR.bed

### Get fasta file
bedtools getfasta -name -s -split -fi reference/fasta/hg38.fa \
-bed /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_3UTR.bed \
> /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_3UTR.fa

bedtools getfasta -name -s -split -fi reference/fasta/hg38.fa \
-bed /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_non-3UTR.bed \
> /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_non-3UTR.fa


# Get region infor for each gene (5UTR, ORF, 3UTR)
python tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/F_split_into_each_section.py \
<(grep -v -w "NA" /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.bed) \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.bed

grep -w "NA" /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.bed \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" "exon" "\t" $4}' \
| LANG=C sort -V -k 1,1 -k 2,2n -k 6,6 >  /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.noncoding.bed

python tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/F_split_into_each_section.py \
<(grep -w "NA" /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $2 "\t" $3 "\t" $9 "\t" $10 "\t" $11 "\t" $12}') \
GRCh38_gencode38_classification.filtered_lite_each_region.noncoding.bed
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "exon" "\t" $6}' GRCh38_gencode38_classification.filtered_lite_each_region.noncoding.bed > /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.noncoding.bed

# Add intron region
python tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/G_add_intron_infor.py  \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.bed \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region_intron_plus.bed

python tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/G_add_intron_infor.py  \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.noncoding.bed \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region_intron_plus.noncoding.bed

### Creating a Markov Background model
## http://meme-suite.org/doc/fasta-get-markov.html
fasta-get-markov -dna \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_3UTR.fa \
/path/to/RipSeq/GSM944520/MEME/markov_background_for_MEME_3UTR.txt

fasta-get-markov -dna \
/path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_non-3UTR.fa \
/path/to/RipSeq/GSM944520/MEME/markov_background_for_MEME_non-3UTR.txt

/path/to/macs3 callpeak -f BAM -t /path/to/RipSeq/${GEO_id}/STAR/merged_STAR_PB29_genome_rm_repbase_rrna/genome_uniq.bam -g hs -B -q 0.05 --nomodel --nolambda --keep-dup all --call-summits --outdir /path/to/RipSeq/${GEO_id}/macs -n genome_uniq

cat /path/to/RipSeq/GSM944520/macs/genome_uniq_peaks.narrowPeak \
| grep "^chr" | grep -v "^chrM"
# 5th: integer score for display. It's calculated as int(-10*log10pvalue) or int(-10*log10qvalue) depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff. Please note that currently this value might be out of the [0-1000] range defined in UCSC ENCODE narrowPeak format. You can let the value saturated at 1000 (i.e. p/q-value = 10^-100) by using the following 1-liner awk: awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak
# 7th: fold-change at peak summit
# 8th: -log10pvalue at peak summit
# 9th: -log10qvalue at peak summit
# 10th: relative summit position to peak start

# → extract high confident peaks in 3'UTR of any isoform (allowing overlapping other 5'UTR/CDS/Intron of other isoforms) using Rscript （→/path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_3UTR.bed)
wc -l /path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_3UTR.bed # 7935 

### make a fasta file for MEME
bedtools getfasta -name -s -split -fi reference/fasta/hg38.fa \
-bed /path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_3UTR.bed \
> /path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_3UTR.fa

bedtools getfasta -name -s -split -fi reference/fasta/hg38.fa \
-bed /path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_non-3UTR.bed \
> /path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_non-3UTR.fa

meme \
/path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_3UTR.fa \
-dna -bfile /path/to/RipSeq/GSM944520/MEME/markov_background_for_MEME_3UTR.txt \
-oc /path/to/RipSeq/GSM944520/MEME -minw 6 -maxw 50 -nmotifs 3 -maxsize 1000000000 -mod zoops -revcomp

meme \
/path/to/RipSeq/GSM944520/macs/LCL_ELAVL1_non-3UTR.fa \
-dna -bfile /path/to/RipSeq/GSM944520/MEME/markov_background_for_MEME_non-3UTR.txt \
-oc /path/to/RipSeq/GSM944520/MEME_non3UTR -minw 6 -maxw 50 -nmotifs 3 -maxsize 1000000000 -mod zoops -revcomp

