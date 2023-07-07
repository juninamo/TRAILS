#!/bin/sh
#$ -S /bin/sh

# Remove 3' end with quality score (-t) below 20

# EvoImmunoPop
ls /path/to/EGA/EGAD00001002714/*/*fq.gz > /path/to/EGA/inamo/fastq_sample.path
while read line
do
    echo "trim_galore ${line}"
    time trim_galore --fastqc -o /path/to/EGA/inamo/trim_galore/ ${line}
done

# ImmVar
cat /path/to/ImmVar/SRA_IID_table.txt | cut -f 1 | while read line; do
    echo "trim_galore ${line}"
    time trim_galore -q 20 -o /path/to/ImmVar/trim_galore/ --paired /path/to/ImmVar/dump/${line}/${line}_dbGaP-24937_1.fastq.gz /path/to/ImmVar/dump/${line}/${line}_dbGaP-24937_2.fastq.gz
done

# DICE
cat /path/to/DICE/cell_type.list | while read cell; do
cat /path/to/DICE/${cell}_sra.list | while read sra; do
    echo "trim_galore ${cell}; ${sra}"
    time trim_galore -q 20 -o /path/to/DICE/trim_galore/${cell}/ /path/to/DICE/dump/${cell}/${sra}/${sra}_dbGaP-24937.fastq.gz
done
done
