#!/bin/bash
#$ -S /bin/sh

# RBPmap

# Get region infor for each gene (5UTR, ORF, 3UTR)
## 1つ目の引数: 各遺伝子の代表Isoformを抽出したBEDファイルを指定。
## 2つ目の引数: 5' UTR、CDS、3' UTRごとに分割したBEDファイルを出力。
python tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/F_split_into_each_section.py \
<(grep -v -w "NA" /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.bed) \
/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.bed

# Add intron region
python tools/NGS-Tutorial/PAR-CLIP_Tutorial_CstF64/G_add_intron_infor.py  \
/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region.bed \
/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region_intron_plus.bed

# for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
# /data02/home/juninamo/tools/RBPmap_1.2/UCSC/faToNib chr${chr}.fa chr${chr}.nib
# done

mkdir -p ${HOME}/tools/RBPmap_1.2/data
mkdir -p ${HOME}/tools/RBPmap_1.2/results
chmod u+x ${HOME}/tools/RBPmap_1.2/results

LANG=C sort -V -k 1,1 -k 2,2n /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_each_region_intron_plus.bed | awk '{print $1 ":" $2 "-" $3 ":" $4}' | uniq > ${HOME}/tools/RBPmap_1.2/data/PB29_exons_introns_coordinate.txt
shuf ${HOME}/tools/RBPmap_1.2/data/PB29_exons_introns_coordinate.txt | split -l 4999 -d - ${HOME}/tools/RBPmap_1.2/data/PB29_exons_introns_coordinate_
ls ${HOME}/tools/RBPmap_1.2/data/PB29_exons_introns_coordinate_* | sed -e "s/\/data02\/home\/juninamo\/tools\/RBPmap_1.2\/data\/PB29_exons_introns_coordinate_//g" > file.list

nohup \
cat file.list | while read line; do
${HOME}/usr/local/perl/bin/perl ${HOME}/tools/RBPmap_1.2/RBPmap.pl -input ${HOME}/tools/RBPmap_1.2/data/PB29_exons_introns_coordinate_${line} -genome human -db hg38 -db_motifs all_human -job_name ${HOME}/tools/RBPmap_1.2/results/PB29_exons_introns_${line}
done > ${HOME}/tools/RBPmap_1.2/log/`date "+%Y%m%d"`.log &

