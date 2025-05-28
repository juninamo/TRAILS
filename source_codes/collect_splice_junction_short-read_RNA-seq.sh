#!/bin/sh
#$ -S /bin/sh

## SJ_files/*_${group}_2nd_SJ.out.tab: output of STAR
## remove junction read =< 3 according to each cell-type

for group in CD4 CD8 B Mono NK PB; do
echo -e ${group}
mkdir --parents /path/to/annotation/SJ_files/${group}
ls /path/to/annotation/SJ_files/*_${group}_2nd_SJ.out.tab | wc -l
awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/annotation/SJ_files/*_${group}_2nd_SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 2 {print $0}' \
> /path/to/annotation/SJ_files/${group}/sum_SJ_thresh3.out.tab

awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/annotation/SJ_files/*_${group}_2nd_SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 0 {print $0}' \
> /path/to/annotation/SJ_files/${group}/sum_SJ.out.tab
done

tac /data01/DICE/inamo/cell_type.list | while read dice; do
echo -e ${dice}
mkdir --parents /path/to/annotation/SJ_files/${dice}
ls /path/to/DICE/BAM/${dice}/*/SJ.out.tab | wc -l
awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/DICE/BAM/${dice}/*/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 2 {print $0}' \
> /path/to/annotation/SJ_files/${dice}/sum_SJ_thresh3.out.tab

awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/DICE/BAM/${dice}/*/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 0 {print $0}' \
> /path/to/annotation/SJ_files/${dice}/sum_SJ.out.tab
done

for i in {1..5}; do
echo "Mono-${i}"
mkdir --parents /path/to/annotation/SJ_files/Mono-${i}
ls /path/to/EGA/BAM/*-${i}_STAR_GRCh38/SJ.out.tab | wc -l
awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/EGA/BAM/*-${i}_STAR_GRCh38/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 2 {print $0}' \
> /path/to/annotation/SJ_files/Mono-${i}/sum_SJ_thresh3.out.tab

awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/EGA/BAM/*-${i}_STAR_GRCh38/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 0 {print $0}' \
> /path/to/annotation/SJ_files/Mono-${i}/sum_SJ.out.tab
done

for i in CD4T-Activated MoDC-unstim MoDC-FLU MoDC-IFNb; do
echo "${i}"
mkdir --parents /path/to/annotation/SJ_files/${i}
ls /path/to/ImmVar/BAM/*-${i}_STAR_GRCh38/SJ.out.tab | wc -l
awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/ImmVar/BAM/*-${i}_STAR_GRCh38/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 2 {print $0}' \
> /path/to/annotation/SJ_files/${i}/sum_SJ_thresh3.out.tab

awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/ImmVar/BAM/*-${i}_STAR_GRCh38/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 0 {print $0}' \
> /path/to/annotation/SJ_files/${i}/sum_SJ.out.tab
done

echo "GEUV, CEU"
mkdir --parents /path/to/annotation/SJ_files/GEUV
ls /path/to/GEUV/BAM/*_STAR_GRCh38/SJ.out.tab | wc -l
awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/GEUV/BAM/*_STAR_GRCh38/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 2 {print $0}' \
> /path/to/annotation/SJ_files/GEUV/sum_SJ_thresh3.out.tab
awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
/path/to/GEUV/BAM/*_STAR_GRCh38/SJ.out.tab \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 0 {print $0}' \
> /path/to/annotation/SJ_files/GEUV/sum_SJ.out.tab

cat /path/to/annotation/SJ_files/*/sum_SJ_thresh3.out.tab \
| LANG=C sort -V -k 1,1 -k 2,2n \
| awk '{ if($0 != "") { a[$1 "_" $2 "_" $3 "_" $4 "_" $5 "_" $6]+=$7; } }END{for(i in a)print i "\t" a[i] "\t" $8 "\t" $9;}' \
| tr "_" "\t" | LANG=C sort -V -k 1,1 -k 2,2n | awk '$7 > 2 {print $0}' \
> /path/to/annotation/SJ_files/all-sr_sum_SJ_thresh3.out.tab


# split by chrosomes
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
echo "chr${chr}"
fgrep -w "chr${chr}" /path/to/annotation/SJ_files/all-sr_sum_SJ_thresh3.out.tab > /path/to/annotation/SJ_files/all-sr_sum_SJ_thresh3_chr${chr}.out.tab
done

