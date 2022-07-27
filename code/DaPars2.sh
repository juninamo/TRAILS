#!/bin/bash
##$ -S /bin/sh

cat /path/to/GEUV/EUR373_sample.txt | cut -f 1 | while read sample; do
echo "bedgraph ${sample}"
bamtocov --wig 0 /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam > /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}.wig
done

cat /path/to/GEUV/YRI89_sample.txt | cut -f 1 | while read sample; do
echo "bedgraph ${sample}"
bamtocov --wig 0 /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam > /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}.wig
done

cat /path/to/DICE/cell_type.list | while read cell; do
cat /path/to/DICE/${cell}_sra.list | while read sra; do
if [ ! -e "/path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}.wig" ]; then
echo "bedgraph ${cell}; ${sra}"
bamtocov --wig 0 /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam > /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}.wig
fi
done
done

cat /path/to/EGA/fq_sample.list | while read line; do
echo "bedgraph ${line}"
bamtocov --wig 0 /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam > /path/to/EGA/BAM/${line}_STAR_PB29/${line}.wig
done

cat /path/to/ImmVar/SRA_IID_table.txt | cut -f 1 | while read line; do
LINE=$(cat /path/to/ImmVar/SRA_IID_table.txt | grep -w ${line} | awk '{print $2}')
echo "bedgraph ${line}; ${LINE}"
bamtocov --wig 0 /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam > /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}.wig
done

mkdir tools/DaPars2/data/wig
ln -s /path/to/*/BAM/*/*_STAR_PB29/*.wig tools/DaPars2/data/wig/
ln -s /path/to/*/BAM/*_STAR_PB29/*.wig tools/DaPars2/data/wig/

# make "/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.coding.bed"
## /path/to/scripts/nanopore/dapars.R
### 1. cluster isoforms by identical coding-sequence using VSEARCH
### 2. make bed file including isoforms predicted coding by all SQANTI3,CPAT and RNAsamba 

python tools/DaPars2/src/DaPars_Extract_Anno_fiveUTR.py \
-b /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.coding.bed \
-s tools/DaPars2/data/read2gene.txt \
-o tools/DaPars2/data/LR_extracted_5UTR.bed

python tools/DaPars2/src/DaPars_Extract_Anno.py \
-b /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.coding.bed \
-s tools/DaPars2/data/read2gene.txt \
-o tools/DaPars2/data/LR_extracted_3UTR.bed


mkdir -p tools/DaPars2/data/result

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo -e "chr${chr}"
cd /path/to/tools/DaPars2/data
python ../src/Dapars2_Multi_Sample_5UTR.py \
config_5UTR.txt \
chr${chr}
done

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo -e "chr${chr}"
cd /path/to/tools/DaPars2/data
python ../src/Dapars2_Multi_Sample.py \
config_3UTR.txt \
chr${chr}
done


