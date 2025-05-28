#!/bin/bash
#$ -S /bin/sh

# Calderon ATAC-seq
## Nat Genet 51, 1494–1505 (2019).

mkdir --parents /path/to/ATAC_seq/Calderon/log
#QC
mkdir --parents /path/to/ATAC_seq/Calderon/trim_galore
while read line
do
    echo "trim_galore ${line}"
  time trim_galore -q 20 --fastqc -o /path/to/ATAC_seq/Calderon/trim_galore/ --paired /path/to/ATAC_seq/Calderon/dump/${line}/${line}_1.fastq.gz /path/to/ATAC_seq/Calderon/dump/${line}/${line}_2.fastq.gz
done < /path/to/ATAC_seq/Calderon/sra_file.list > /path/to/ATAC_seq/Calderon/log/`date "+%Y%m%d"`_trim_galore.log 2>&1

tools/bowtie2-2.2.9/bowtie2-build -f reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna tools/bowtie2-2.2.9/index/GRCh38

cat /path/to/ATAC_seq/Calderon/sra_file.list | while read line; do
echo "bowtie2 ${line}"
tools/bowtie2-2.2.9/bowtie2 -p 30 -x tools/bowtie2-2.2.9/index/GRCh38 \
-1 /path/to/ATAC_seq/Calderon/trim_galore/${line}_1_val_1.fq.gz \
-2 /path/to/ATAC_seq/Calderon/trim_galore/${line}_2_val_2.fq.gz \
-S /path/to/ATAC_seq/Calderon/BAM/${line}_GRCh38.sam
done > /path/to/ATAC_seq/Calderon/log/`date "+%Y%m%d"`_bowtie_GRCh38.log 2>&1

cat /path/to/ATAC_seq/Calderon/sra_file.list | while read line; do
    if [ ! -e "/path/to/ATAC_seq/Calderon/BAM/${line}_sorted_GRCh38.bam.bai" ]; then
    echo "samtools ${line}"
    /path/to/samtools sort -@ 10 -O bam -o /path/to/ATAC_seq/Calderon/BAM/${line}_sorted_GRCh38.bam /path/to/ATAC_seq/Calderon/BAM/${line}_GRCh38.sam # sam を sort しながら bam へ変換
    /path/to/samtools index -@ 10 /path/to/ATAC_seq/Calderon/BAM/${line}_sorted_GRCh38.bam
    fi
done

cat /path/to/ATAC_seq/Calderon/sra_file.list | while read line; do
    echo "samtools ${line}"
    /path/to/samtools view -b -q 10 /path/to/ATAC_seq/Calderon/BAM/${line}_sorted_GRCh38.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > /path/to/ATAC_seq/Calderon/BAM/${line}_sorted_filtered_GRCh38.bam
    java -jar tools/picard-2.25.0/picard.jar MarkDuplicates I=/path/to/ATAC_seq/Calderon/BAM/${line}_sorted_filtered_GRCh38.bam O=/path/to/ATAC_seq/Calderon/BAM/${line}_sorted_remDup_GRCh38.bam M=/path/to/ATAC_seq/Calderon/BAM/${line}_dups_GRCh38.txt REMOVE_DUPLICATES=true
    /path/to/samtools index -@ 10 /path/to/ATAC_seq/Calderon/BAM/${line}_sorted_remDup_GRCh38.bam
done

cat /path/to/ATAC_seq/Calderon/sra_file.list | while read line; do
    LINE=$(cat /path/to/ATAC_seq/Calderon/sample_summary.txt | grep ${line} | awk '{print $4}')
    echo $LINE: $line
    macs3 callpeak -f BAMPE -t /path/to/ATAC_seq/Calderon/BAM/${line}_sorted_remDup_GRCh38.bam -g hs -B -q 0.05 --nomodel --nolambda --keep-dup all --call-summits --outdir /path/to/ATAC_seq/Calderon/macs/GRCh38 -n ${line}_${LINE}
done > /path/to/ATAC_seq/Calderon/log/`date "+%Y%m%d"`_macs3.log 2>&1

# merge by cell-types
sed -e "1d" /path/to/ATAC_seq/Calderon/sample_summary.txt | cut -f 4 | sort | uniq > cell_types.list
cat cell_types.list | while read cell; do
    # no_treatment
    mkdir --parents tmp/${cell}
    echo -e "no_treatment: copy all ${cell} files in tmp/${cell}..."
    grep -w ${cell} /path/to/ATAC_seq/Calderon/sample_summary.txt | grep "no_treament" | awk '{print $1}' > sample.list
    cat sample.list | while read sample; do
        cp /path/to/ATAC_seq/Calderon/BAM/${sample}_sorted_remDup_GRCh38.bam tmp/${cell}/
    done
    echo -e "no_treatment: merging all ${cell} files..."
    /path/to/samtools merge -@ 10 tmp/${cell}/finalBamFile.bam tmp/${cell}/*.bam
    echo -e "no_treatment: sorting ${cell} file..."
    /path/to/samtools sort -@ 10 tmp/${cell}/finalBamFile.bam -o tmp/${cell}/${cell}_no_treament_merged.sorted.bam
    echo -e "no_treatment: indexing ${cell} file..."
    /path/to/samtools index -@ 10 tmp/${cell}/${cell}_no_treament_merged.sorted.bam
    cp tmp/${cell}/${cell}_no_treament_merged.sorted.bam* /path/to/ATAC_seq/Calderon/BAM/
    echo -e "no_treatment: callpeak ${cell} file..."
    /path/to/macs3 callpeak -f BAMPE -t tmp/${cell}/${cell}_no_treament_merged.sorted.bam -g hs -B -q 0.05 --nomodel --nolambda --keep-dup all --call-summits --outdir /path/to/ATAC_seq/Calderon/macs/GRCh38 -n ${cell}_no_treament_merged
    rm -r tmp/${cell}
    BAMscale scale --bam /path/to/ATAC_seq/Calderon/BAM/${cell}_no_treament_merged.sorted.bam -t 30 --outdir /path/to/ATAC_seq/Calderon/macs/GRCh38
    
    # treatment
    mkdir --parents tmp/${cell}
    echo -e "treatment: copy all ${cell} files in tmp/${cell}..."
    grep -w ${cell} /path/to/ATAC_seq/Calderon/sample_summary.txt | grep -v "no_treament" | awk '{print $1}' > sample.list
    cat sample.list | while read sample; do
        cp /path/to/ATAC_seq/Calderon/BAM/${sample}_sorted_remDup_GRCh38.bam tmp/${cell}/
    done
    echo -e "treatment: merging all ${cell} files..."
    /path/to/samtools merge -@ 10 tmp/${cell}/finalBamFile.bam tmp/${cell}/*.bam
    echo -e "treatment: sorting ${cell} file..."
    /path/to/samtools sort -@ 10 tmp/${cell}/finalBamFile.bam -o tmp/${cell}/${cell}_treament_merged.sorted.bam
    echo -e "treatment: indexing ${cell} file..."
    /path/to/samtools index -@ 10 tmp/${cell}/${cell}_treament_merged.sorted.bam
    cp tmp/${cell}/${cell}_treament_merged.sorted.bam* /path/to/ATAC_seq/Calderon/BAM/
    echo -e "treatment: callpeak ${cell} file..."
    /path/to/macs3 callpeak -f BAMPE -t tmp/${cell}/${cell}_treament_merged.sorted.bam -g hs -B -q 0.05 --nomodel --nolambda --keep-dup all --call-summits --outdir /path/to/ATAC_seq/Calderon/macs/GRCh38 -n ${cell}_treament_merged
    rm -r tmp/${cell}
    BAMscale scale --bam /path/to/ATAC_seq/Calderon/BAM/${cell}_treament_merged.sorted.bam -t 30 --outdir /path/to/ATAC_seq/Calderon/macs/GRCh38
done

# tabix and index
cat cell_types.list | while read cell; do
echo "${cell}"
bgzip -@ 8 /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_treament_merged_treat_pileup.bdg
tabix -p bed /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_treament_merged_treat_pileup.bdg.gz
bgzip -@ 8 /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_no_treament_merged_treat_pileup.bdg
tabix -p bed /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_no_treament_merged_treat_pileup.bdg.gz
done

cat cell_types.list | while read cell; do
echo "${cell}"
zcat /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_no_treament_merged_treat_pileup.bdg.gz | awk '{sum+=$4;} END{print sum;}' > /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_no_treament_merged_treat_pileup.bdg.gz.sum
zcat /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_treament_merged_treat_pileup.bdg.gz | awk '{sum+=$4;} END{print sum;}' > /path/to/ATAC_seq/Calderon/macs/GRCh38/${cell}_treament_merged_treat_pileup.bdg.gz.sum
done
