#!/bin/sh
#$ -S /bin/sh

# remapping to isoform atlas by STAR

mkdir --parents tools/STAR/genome_GRCh38_PB29
/opt/bin/STAR --runMode genomeGenerate \
--genomeDir tools/STAR/genome_GRCh38_PB29/ \
--genomeFastaFiles reference/GENCODE/GRCh38/release_38/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /path/to/GRCh38_gencode38_classification.filtered_lite.gtf \
--sjdbOverhang 100 \
--runThreadN 35

cat /path/to/GEUV/EUR373_sample.txt | cut -f 1 | while read sample; do
    echo "mapping ${sample}"
    mkdir --parents /path/to/GEUV/BAM/${sample}_STAR_PB29
    /opt/bin/STAR --outFilterType BySJout \
    --runThreadN 10 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir tools/STAR/genome_GRCh38_PB29/ \
    --sjdbGTFfile /path/to/GRCh38_gencode38_classification.filtered_lite.gtf \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix /path/to/GEUV/BAM/${sample}_STAR_PB29/ \
    --readFilesCommand gunzip -c \
    --readFilesIn /opt/mnt/QNAPGFD01_fastq/GEUV/CEU/${sample}_R1_QC.fastq.gz /opt/mnt/QNAPGFD01_fastq/GEUV/CEU/${sample}_R2_QC.fastq.gz

    echo "samtools ${sample}"
    /opt/bin/samtools sort -@ 10 /path/to/GEUV/BAM/${sample}_STAR_PB29/Aligned.sortedByCoord.out.bam -o /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam
    rm /path/to/GEUV/BAM/${sample}_STAR_PB29/Aligned.sortedByCoord.out.bam
    /opt/bin/samtools index -@ 10 /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam
done

cat /path/to/DICE/cell_type.list | while read cell; do
mkdir --parents /path/to/DICE/BAM/${cell}
cat /path/to/DICE/${cell}_sra.list | while read sra; do
mkdir --parents /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29
    echo "mapping ${cell}; ${sra}"
    /opt/bin/STAR --outFilterType BySJout \
    --runThreadN 10 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir tools/STAR/genome_GRCh38_PB29/ \
    --sjdbGTFfile /path/to/GRCh38_gencode38_classification.filtered_lite.gtf \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/ \
    --readFilesCommand gunzip -c \
    --readFilesIn /path/to/DICE/trim_galore/${cell}/${sra}_dbGaP-24937_trimmed.fq.gz

    echo "samtools ${cell}; ${sra}"
    /opt/bin/samtools sort -@ 10 /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/Aligned.sortedByCoord.out.bam -o /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam
    rm /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/Aligned.sortedByCoord.out.bam
    /opt/bin/samtools index -@ 10 /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam
done
done

cat /path/to/EGA/fq_sample.list | while read line; do
    echo "mapping ${line}"
    mkdir --parents /path/to/EGA/BAM/${line}_STAR_PB29
    /opt/bin/STAR --outFilterType BySJout \
    --runThreadN 10 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir tools/STAR/genome_GRCh38_PB29/ \
    --sjdbGTFfile /path/to/GRCh38_gencode38_classification.filtered_lite.gtf \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix /path/to/EGA/BAM/${line}_STAR_PB29/ \
    --readFilesCommand gunzip -c \
    --readFilesIn /path/to/EGA/trim_galore/${line}_aligned_trimmed.fq.gz

    echo "samtools ${line}"
    /opt/bin/samtools sort -@ 10 /path/to/EGA/BAM/${line}_STAR_PB29/Aligned.sortedByCoord.out.bam -o /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam
    rm /path/to/EGA/BAM/${line}_STAR_PB29/Aligned.sortedByCoord.out.bam
    /opt/bin/samtools index -@ 10 /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam
done

cat /path/to/ImmVar/SRA_IID_table.txt | cut -f 1 | while read line; do
LINE=$(cat /path/to/ImmVar/SRA_IID_table.txt | grep -w ${line} | awk '{print $2}')
    echo "mapping ${line}; ${LINE}"
    mkdir --parents /path/to/ImmVar/BAM/${LINE}_STAR_PB29
    /opt/bin/STAR --outFilterType BySJout \
    --runThreadN 10 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir tools/STAR/genome_GRCh38_PB29/ \
    --sjdbGTFfile /path/to/GRCh38_gencode38_classification.filtered_lite.gtf \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix /path/to/ImmVar/BAM/${LINE}_STAR_PB29/ \
    --readFilesCommand gunzip -c \
    --readFilesIn /path/to/ImmVar/trim_galore/${line}_dbGaP-24937_1_val_1.fq.gz /path/to/ImmVar/trim_galore/${line}_dbGaP-24937_2_val_2.fq.gz

    echo "samtools ${line}; ${LINE}"
    /opt/bin/samtools sort -@ 10 /path/to/ImmVar/BAM/${LINE}_STAR_PB29/Aligned.sortedByCoord.out.bam -o /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam
    rm /path/to/ImmVar/BAM/${LINE}_STAR_PB29/Aligned.sortedByCoord.out.bam
    /opt/bin/samtools index -@ 10 /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam
done