#!/bin/sh
#$ -S /bin/sh

# mapping RNA-seq (hapmap)
for rna_id in SRR031811 SRR031815 SRR031820 SRR031821 SRR031826 SRR031827 SRR031829 SRR031832 SRR031834 SRR031836 SRR031840 SRR031843 SRR031846 SRR031847 SRR031848 SRR031849 SRR031852 SRR031854 SRR031856 SRR031857 SRR031858 SRR031859 SRR031861 SRR031862 SRR031863 SRR031865 SRR031867 SRR031869 SRR031873 SRR031874 SRR031875 SRR031876 SRR031878 SRR031880 SRR031882 SRR031884 SRR031887 SRR031889 SRR031891 SRR031893 SRR031894 SRR031896 SRR031897 SRR031900 SRR031902 SRR031903 SRR031905 SRR031908 SRR031910 SRR031911 SRR031915 SRR031916 SRR031917 SRR031921 SRR031923 SRR031924 SRR031925 SRR031927 SRR031929 SRR031930 SRR031933 SRR031935 SRR031936 SRR031938 SRR031939 SRR031942 SRR031945 SRR031947 SRR031950 SRR031951 SRR031953 SRR031955 SRR031958 SRR031959 SRR031961 SRR031963 SRR031968 SRR031969 SRR031971 SRR031812 SRR031813 SRR031814 SRR031816 SRR031817 SRR031818 SRR031819 SRR031822 SRR031823 SRR031824 SRR031825 SRR031828 SRR031830 SRR031831 SRR031833 SRR031835 SRR031837 SRR031838 SRR031839 SRR031841 SRR031842 SRR031844 SRR031845 SRR031850 SRR031851 SRR031853 SRR031855 SRR031860 SRR031864 SRR031866 SRR031868 SRR031870 SRR031871 SRR031872 SRR031877 SRR031879 SRR031881 SRR031883 SRR031885 SRR031886 SRR031888 SRR031890 SRR031892 SRR031895 SRR031898 SRR031899 SRR031901 SRR031904 SRR031906 SRR031907 SRR031909 SRR031912 SRR031913 SRR031914 SRR031918 SRR031919 SRR031920 SRR031922 SRR031926 SRR031928 SRR031931 SRR031932 SRR031934 SRR031937 SRR031940 SRR031941 SRR031943 SRR031944 SRR031946 SRR031948 SRR031949 SRR031952 SRR031954 SRR031956 SRR031957 SRR031960 SRR031962 SRR031964 SRR031965 SRR031966 SRR031967 SRR031970; do
rna_bam="/path/to/hapmap/RNAseq/GSE19480/STAR/${rna_id}_STAR_PB29/uniq_sorted.bam"
stringtie_out="/path/to/hapmap/RNAseq/GSE19480/StringTie_PB29/${rna_id}"
transcript_gtf="/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf"
mkdir --parents ${stringtie_out}
stringtie -e -B -p 10 \
-G ${transcript_gtf} \
-o ${stringtie_out}/${rna_id}_stringtie.gtf \
${rna_bam}
done

# mapping Ribo-seq (hapmap)
for ribo_id in SRR1585486 SRR1585487 SRR1585488 SRR1585490 SRR1585492 SRR1585494 SRR1585495 SRR1585496 SRR1585498 SRR1585499 SRR1585500 SRR1585507 SRR1585508 SRR1585510 SRR1585511 SRR1585512 SRR1585513 SRR1585515 SRR1585516 SRR1585517 SRR1585518 SRR1585519 SRR1585521 SRR1585522 SRR1585523 SRR1585524 SRR1585525 SRR1585528 SRR1585529 SRR1585530 SRR1585531 SRR1585533 SRR1585534 SRR1585535 SRR1585536 SRR1585537 SRR1585538 SRR1585539 SRR1585540 SRR1585541 SRR1585542 SRR1585543 SRR1585546 SRR1585547 SRR1585548 SRR1585549 SRR1585550 SRR1585551 SRR1585552 SRR1585553 SRR1585554 SRR1585557; do
trim_out="/path/to/riboseq/trim_galore/GSE61742/"
riboseq_fq="/path/to/riboseq/GSE61742/dump/${ribo_id}/${ribo_id}.fastq.gz"
rrna_idx="/path/to/tools/STAR/rRNA_contam"
transcript_fa="/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fasta"
transcript_idx="/path/to/tools/STAR/genome_GRCh38_PB29_transcript"
transcript_gtf="/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf"
star_ribo_out="/path/to/riboseq/STAR/GSE61742/"
stringtie_out="/path/to/riboseq/StringTie_PB29/GSE61742/${ribo_id}"
ribo_bam="${star_ribo_out}${ribo_id}_transcriptome_Aligned.sortedByCoord.out.bam"
nproc=8 # threads

# Ribo-seq: to remove rRNA
## https://imamachi-n.hatenablog.com/entry/2017/03/31/215218
mkdir -p tools/STAR/rRNA_contam
/path/to/STAR \
--runThreadN ${nproc} \
--runMode genomeGenerate \
--genomeDir ${rrna_idx} \
--genomeFastaFiles reference/fasta/contam_Ribosomal_RNA.fa

# Ribo-seq: mapping to transcriptome
mkdir -p ${transcript_idx}
/path/to/STAR \
--runThreadN ${nproc} \
--runMode genomeGenerate \
--genomeDir ${transcript_idx} \
--genomeSAindexNbases 11 --genomeChrBinNbits 12 \
--genomeFastaFiles ${transcript_fa}


#============================================
# step 0: trim adaptor sequence
#============================================

trim_galore -q 20 -o ${trim_out} ${riboseq_fq}

#============================================
# step 1: filter rrna
#============================================
/path/to/STAR \
--runMode alignReads \
--runThreadN $nproc \
--genomeDir $rrna_idx \
--readFilesIn ${trim_out}${ribo_id}_trimmed.fq.gz \
--readFilesCommand gunzip -c \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix ${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastq \
--outSAMattributes All \
--outStd BAM_Unsorted \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--alignEndsType EndToEnd > ${star_ribo_out}${ribo_id}_repbase_rrna_comtam.bam

#============================================
# step 2: mapping to transcriptome
#============================================
/path/to/STAR \
--runMode alignReads \
--runThreadN $nproc \
--genomeDir tools/STAR/genome_GRCh38_PB29_transcript \
--readFilesIn ${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastq \
--outFilterMultimapNmax 8 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
--outFilterMismatchNmax 4 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFileNamePrefix ${star_ribo_out}${ribo_id}_transcriptome_ \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 32000000000 \
--outFilterType BySJout \
--outReadsUnmapped Fastx

/path/to/samtools sort -@ 8 ${star_ribo_out}${ribo_id}_transcriptome_Aligned.sortedByCoord.out.bam -o ${star_ribo_out}${ribo_id}_transcriptome_sorted.bam
/path/to/samtools index -@ 8 ${star_ribo_out}${ribo_id}_transcriptome_sorted.bam
rm ${star_ribo_out}${ribo_id}_transcriptome_Aligned.sortedByCoord.out.bam

mkdir --parents ${stringtie_out}
stringtie -e -B -p 10 \
-G ${transcript_gtf} \
-o ${stringtie_out}/${ribo_id}_stringtie.gtf \
${star_ribo_out}${ribo_id}_transcriptome_sorted.bam

done


# mapping RNA-seq (GEUVADIS)
cat /path/to/GEUV/sample.txt | cut -f 1 | while read sample; do
/path/to/STAR --outFilterType BySJout \
--runThreadN 10 \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--genomeDir tools/STAR/genome_GRCh38_PB29/ \
--sjdbGTFfile /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
--outSAMstrandField intronMotif \
--outFileNamePrefix /path/to/GEUV/BAM/${sample}_STAR_PB29/ \
--readFilesCommand gunzip -c \
--readFilesIn /path/to/GEUV/CEU/${sample}_R1_QC.fastq.gz /path/to/GEUV/CEU/${sample}_R2_QC.fastq.gz
echo "samtools ${sample}"
/path/to/samtools sort -@ 10 /path/to/GEUV/BAM/${sample}_STAR_PB29/Aligned.sortedByCoord.out.bam -o /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam
rm /path/to/GEUV/BAM/${sample}_STAR_PB29/Aligned.sortedByCoord.out.bam
/path/to/samtools index -@ 10 /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam
done
