#!/bin/sh
#$ -S /bin/sh

# by kallisto

kallisto index -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered --make-unique /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fa

# lengh
## DICE → 50 bp
## EGA → 100 bp

## ImmVar
cat /path/to/ImmVar/SRA_IID_table.txt | cut -f 1 | while read line; do
    
    LINE=$(cat /path/to/ImmVar/SRA_IID_table.txt | grep -w ${line} | awk '{print $2}')
    if [ ! -e "/path/to/ImmVar/kallisto/${LINE}/abundance.tsv" ]; then

    echo -e "kallisto: ${LINE} ..."
    mkdir --parents /path/to/ImmVar/kallisto/${LINE}
    result_dir="/path/to/ImmVar/kallisto/${LINE}"
    kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered -o ${result_dir} /path/to/ImmVar/trim_galore/${line}_dbGaP-24937_1_val_1.fq.gz /path/to/ImmVar/trim_galore/${line}_dbGaP-24937_2_val_2.fq.gz

    fi
done

## GEUV
cat /path/to/GEUV/EUR373_sample.txt | cut -f 1 | while read sample; do
    if [ ! -e "/path/to/GEUV/kallisto/${sample}/abundance.tsv" ]; then

    echo -e "kallisto: ${sample} ..."
    mkdir --parents /path/to/GEUV/kallisto/${sample}
    result_dir="/path/to/GEUV/kallisto/${sample}"
    kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered -o ${result_dir} /opt/mnt/QNAPGFD01_fastq/GEUV/CEU/${sample}_R1_QC.fastq.gz /opt/mnt/QNAPGFD01_fastq/GEUV/CEU/${sample}_R2_QC.fastq.gz
    
    fi
done

## EGA
cat /path/to/EGA/fq_sample.list | while read line; do
    if [ ! -e "/path/to/EGA/kallisto/${line}/abundance.tsv" ]; then

    echo -e "kallisto: ${line} ..."
    mkdir --parents /path/to/EGA/kallisto/${line}
    result_dir="/path/to/EGA/kallisto/${line}"
    kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered -o ${result_dir} --single -l 100 -s 20 /path/to/EGA/trim_galore/${line}_aligned_trimmed.fq.gz

    fi
done


## DICE
cat /path/to/DICE/cell_type.list | while read cell; do
    cat /path/to/DICE/${cell}_sra.list | while read sra; do
        if [ ! -e "/path/to/DICE/kallisto/${sra}/abundance.tsv" ]; then

        echo -e "kallisto: ${cell}, ${sra} ..."
        mkdir --parents /path/to/DICE/kallisto/${sra}
        result_dir="/path/to/DICE/kallisto/${sra}"
        kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered -o ${result_dir} --single -l 50 -s 10 /path/to/DICE/trim_galore/${cell}/${sra}_dbGaP-24937_trimmed.fq.gz

        fi
    done
done




## clustered isoforms by VSEARCH

seqkit seq -n /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.vsearch.centroids.fasta | tr "_" "\t" | awk '{print $1}' > id.txt
seqkit grep -f id.txt /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fa -o /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.clustered-by-cds.vsearch.fa

kallisto index -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered_clustered-by-cds-vsearch --make-unique /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.clustered-by-cds.vsearch.fa

## ImmVar
cat /path/to/ImmVar/SRA_IID_table.txt | cut -f 1 | while read line; do
    
    LINE=$(cat /path/to/ImmVar/SRA_IID_table.txt | grep -w ${line} | awk '{print $2}')
    if [ ! -e "/path/to/ImmVar/kallisto_clustered-by-cds-vsearch/${LINE}/abundance.tsv" ]; then

    echo -e "kallisto: ${LINE} ..."
    mkdir --parents /path/to/ImmVar/kallisto_clustered-by-cds-vsearch/${LINE}
    result_dir="/path/to/ImmVar/kallisto_clustered-by-cds-vsearch/${LINE}"
    kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered_clustered-by-cds-vsearch -o ${result_dir} /path/to/ImmVar/trim_galore/${line}_dbGaP-24937_1_val_1.fq.gz /path/to/ImmVar/trim_galore/${line}_dbGaP-24937_2_val_2.fq.gz

    fi
done


## GEUV
cat /path/to/GEUV/EUR373_sample.txt | cut -f 1 | while read sample; do
    if [ ! -e "/path/to/GEUV/kallisto_clustered-by-cds-vsearch/${sample}/abundance.tsv" ]; then

    echo -e "kallisto: ${sample} ..."
    mkdir --parents /path/to/GEUV/kallisto_clustered-by-cds-vsearch/${sample}
    result_dir="/path/to/GEUV/kallisto_clustered-by-cds-vsearch/${sample}"
    kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered_clustered-by-cds-vsearch -o ${result_dir} /opt/mnt/QNAPGFD01_fastq/GEUV/CEU/${sample}_R1_QC.fastq.gz /opt/mnt/QNAPGFD01_fastq/GEUV/CEU/${sample}_R2_QC.fastq.gz
    
    fi
done

## EGA
cat /path/to/EGA/fq_sample.list | while read line; do
    if [ ! -e "/path/to/EGA/kallisto_clustered-by-cds-vsearch/${line}/abundance.tsv" ]; then

    echo -e "kallisto: ${line} ..."
    mkdir --parents /path/to/EGA/kallisto_clustered-by-cds-vsearch/${line}
    result_dir="/path/to/EGA/kallisto_clustered-by-cds-vsearch/${line}"
    kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered_clustered-by-cds-vsearch -o ${result_dir} --single -l 100 -s 20 /path/to/EGA/trim_galore/${line}_aligned_trimmed.fq.gz

    fi
done


## DICE
cat /path/to/DICE/cell_type.list | while read cell; do
    cat /path/to/DICE/${cell}_sra.list | while read sra; do
        if [ ! -e "/path/to/DICE/kallisto_clustered-by-cds-vsearch/${sra}/abundance.tsv" ]; then

        echo -e "kallisto: ${cell}, ${sra} ..."
        mkdir --parents /path/to/DICE/kallisto_clustered-by-cds-vsearch/${sra}
        result_dir="/path/to/DICE/kallisto_clustered-by-cds-vsearch/${sra}"
        kallisto quant -t 35 -i /path/to/kallisto_index_GRCh38_all-sr-correct-filtered_clustered-by-cds-vsearch -o ${result_dir} --single -l 50 -s 10 /path/to/DICE/trim_galore/${cell}/${sra}_dbGaP-24937_trimmed.fq.gz

        fi
    done
done




# by StringTie2

### ImmVar
mkdir --parents /path/to/ImmVar/StringTie_PB29/ballgown
for i in CD4T_Activated MoDC_unstim MoDC_FLU MoDC_IFNb; do
    cat /path/to/ImmVar/${i}_RNAseq_IID.list | while read line; do

        cell=$(sed -e "s/_/-/g" <(echo -e "${i}"))
        LINE=$(cat /path/to/ImmVar/SRA_IID_table.txt | grep -w ${line} | grep ${cell} | awk '{print $2}')

        if [ ! -e "/path/to/ImmVar/StringTie_PB29/ballgown/${LINE}/${LINE}_stringtie.gtf" ] && [ -e "/path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam" ] ; then
        echo -e "StringTie: ${LINE}_stringtie.gtf ..."
        mkdir --parents /path/to/ImmVar/StringTie_PB29/ballgown/${LINE}

        stringtie -e -B -p 10 \
        -G /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
        -o /path/to/ImmVar/StringTie_PB29/ballgown/${LINE}/${LINE}_stringtie.gtf \
        /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam

        fi
    done
done
ls /path/to/ImmVar/StringTie_PB29/ballgown/*/*-*.gtf | awk -F "/" '{print $8}' | sed -e "s/_stringtie.gtf//g" > ImmVar_sample.list
first=$(head -n 1 ImmVar_sample.list)
grep $'StringTie\ttranscript' /path/to/ImmVar/StringTie_PB29/ballgown/${first}/${first}_stringtie.gtf | awk -F ";" '{print $1 "\t" $2 "\t" $3}' | sed -e 's/gene_id "//g' | sed -e 's/transcript_id//g' | sed -e 's/ref_gene_name//g' | sed -e 's/"/\t/g'  | sed -E 's/[\t ]+/\t/g' |  awk '{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $10 "\t" $9}' | sed -e "s/^chr//g" | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > ImmVar_sample_0.txt
grep -v "^chr" ImmVar_sample_0.txt | LANG=C sort | uniq | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > /path/to/ImmVar/StringTie_PB29/geneInfo.txt

### EGA
mkdir --parents /path/to/EGA/StringTie_PB29/ballgown
for i in {1..5}; do
    cat /path/to/EGA/fq_sample.list | grep -e "-${i}$" | while read line; do

        if [ ! -e "/path/to/EGA/StringTie_PB29/ballgown/${line}/${line}_stringtie.gtf" ] && [ -e "/path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam" ] ; then
        echo -e "StringTie: ${line}_stringtie.gtf ..."
        mkdir --parents /path/to/EGA/StringTie_PB29/ballgown/${line}

        stringtie -e -B -p 10 \
        -G /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
        -o /path/to/EGA/StringTie_PB29/ballgown/${line}/${line}_stringtie.gtf \
        /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam

        fi
    done
done
ls /path/to/EGA/StringTie_PB29/ballgown/*/*-*.gtf | awk -F "/" '{print $8}' | sed -e "s/_stringtie.gtf//g" > EGA_sample.list
first=$(head -n 1 EGA_sample.list)
grep $'StringTie\ttranscript' /path/to/EGA/StringTie_PB29/ballgown/${first}/${first}_stringtie.gtf | awk -F ";" '{print $1 "\t" $2 "\t" $3}' | sed -e 's/gene_id "//g' | sed -e 's/transcript_id//g' | sed -e 's/ref_gene_name//g' | sed -e 's/"/\t/g'  | sed -E 's/[\t ]+/\t/g' |  awk '{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $10 "\t" $9}' | sed -e "s/^chr//g" | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > EGA_sample_0.txt
grep -v "^chr" EGA_sample_0.txt | LANG=C sort | uniq | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > /path/to/EGA/StringTie_PB29/geneInfo.txt


### DICE
mkdir --parents /path/to/DICE/StringTie_PB29/ballgown
cat /path/to/DICE/cell_type.list | while read cell; do
    cat /path/to/DICE/${cell}_sra.list | while read sra; do

        if [ ! -e "/path/to/DICE/StringTie_PB29/ballgown/${sra}/${sra}_stringtie.gtf" ] && [ -e "/path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam" ] ; then
        echo -e "StringTie: ${cell}, ${sra}_stringtie.gtf ..."
        mkdir --parents /path/to/DICE/StringTie_PB29/ballgown/${sra}

        stringtie -e -B -p 10 \
        -G /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
        -o /path/to/DICE/StringTie_PB29/ballgown/${sra}/${sra}_stringtie.gtf \
        /path/to/DICE/BAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam

        fi

    done
done
ls /path/to/DICE/StringTie_PB29/ballgown/*/SRR*.gtf | awk -F "/" '{print $8}' | sed -e "s/_stringtie.gtf//g" > DICE_sample.list
first=$(head -n 1 DICE_sample.list)
grep $'StringTie\ttranscript' /path/to/DICE/StringTie_PB29/ballgown/${first}/${first}_stringtie.gtf | awk -F ";" '{print $1 "\t" $2 "\t" $3}' | sed -e 's/gene_id "//g' | sed -e 's/transcript_id//g' | sed -e 's/ref_gene_name//g' | sed -e 's/"/\t/g'  | sed -E 's/[\t ]+/\t/g' |  awk '{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $10 "\t" $9}' | sed -e "s/^chr//g" | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > DICE_sample_0.txt
grep -v "^chr" DICE_sample_0.txt | LANG=C sort | uniq | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > /path/to/DICE/StringTie_PB29/geneInfo.txt

### GEUV
mkdir --parents /path/to/GEUV/StringTie_PB29/ballgown
cat /path/to/GEUV/all_sample.txt | cut -f 1 | while read sample; do

    if [ ! -e "/path/to/GEUV/StringTie_PB29/ballgown/${sample}/${sample}_stringtie.gtf" ] && [ -e "/path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam" ] ; then
    echo -e "StringTie: ${sample}_stringtie.gtf ..."
    mkdir --parents /path/to/GEUV/StringTie_PB29/ballgown/${sample}

    stringtie -e -B -p 10 \
    -G /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
    -o /path/to/GEUV/StringTie_PB29/ballgown/${sample}/${sample}_stringtie.gtf \
    /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam

    fi
done
ls /path/to/GEUV/StringTie_PB29/ballgown/*/ERR*.gtf | awk -F "/" '{print $8}' | sed -e "s/_stringtie.gtf//g" > GEUV_sample.list
first=$(head -n 1 GEUV_sample.list)
grep $'StringTie\ttranscript' /path/to/GEUV/StringTie_PB29/ballgown/${first}/${first}_stringtie.gtf | awk -F ";" '{print $1 "\t" $2 "\t" $3}' | sed -e 's/gene_id "//g' | sed -e 's/transcript_id//g' | sed -e 's/ref_gene_name//g' | sed -e 's/"/\t/g'  | sed -E 's/[\t ]+/\t/g' |  awk '{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $10 "\t" $9}' | sed -e "s/^chr//g" | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > GEUV_sample_0.txt
grep -v "^chr" GEUV_sample_0.txt | LANG=C sort | uniq | sed -e "1i chr\tstart\tend\tstrand\ttranscript_id\tgene_id" > /path/to/GEUV/StringTie_PB29/geneInfo.txt
