#!/bin/sh
#$ -S /bin/sh

# leafcutter

## extract exon and make annotation code
grep -w "exon" /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf | awk '{print $1 " " $4 " " $5 " " $7 " " $12}' | sed -e 's/"//g' | sed -e 's/;//g' > /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.exon.txt
~/tools/leafcutter/leafviz/gtf2leafcutter.pl -o tools/leafcutter/annotation_codes/PB29 /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf


# ImmVar
for i in CD4T_Activated MoDC_unstim MoDC_FLU MoDC_IFNb; do
    mkdir --parents /path/to/ImmVar/leafcutter_PB29/${i}
    cat /path/to/ImmVar/${i}_RNAseq_IID.list | while read line; do        
        cell=$(sed -e "s/_/-/g" <(echo -e "${i}"))
        LINE=$(cat /path/to/ImmVar/SRA_IID_table.txt | grep -w ${line} | grep ${cell} | awk '{print $2}')
        echo -e "regtools junctions extract: ${LINE}_sorted.bam.junc ..."
        regtools junctions extract \
        -s 0 -a 8 -m 50 -M 500000 \
        /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam \
        -o /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam.junc
        cat /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted.bam.junc | grep "^chr" | grep -v "^chrM" > /path/to/ImmVar/BAM/${LINE}_STAR_PB29/${LINE}_sorted_filtered.bam.junc
    done
done



# EGA
for i in {1..5}; do
    mkdir --parents /path/to/EGA/leafcutter_PB29/Mono_${i}
    cat /path/to/EGA/fq_sample.list | grep -e "-${i}$" | while read line; do        
        echo -e "regtools junctions extract: ${line}_sorted.bam.junc ..."
        regtools junctions extract \
        -s 0 -a 8 -m 50 -M 500000 \
        /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam \
        -o /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam.junc
        cat /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted.bam.junc | grep "^chr" | grep -v "^chrM" > /path/to/EGA/BAM/${line}_STAR_PB29/${line}_sorted_filtered.bam.junc
    done
done


# DICE
cat /path/to/DICE/cell_type.list | while read cell; do
    mkdir --parents /path/to/DICEleafcutter_PB29/${cell}
    cat /path/to/DICE/${cell}_sra.list | while read sra; do      
        echo -e "regtools junctions extract: ${cell}, ${sra}_sorted.bam.junc ..."
        regtools junctions extract \
        -s 0 -a 8 -m 50 -M 500000 \
        /path/to/DICEBAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam \
        -o /path/to/DICEBAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam.junc

        cat /path/to/DICEBAM/${cell}/${sra}_STAR_PB29/${sra}_sorted.bam.junc | grep "^chr" | grep -v "^chrM" > /path/to/DICEBAM/${cell}/${sra}_STAR_PB29/${sra}_sorted_filtered.bam.junc
    done
done


# GEUV

cat /path/to/GEUV/all_sample.txt | cut -f 1 | while read sample; do        
    echo -e "regtools junctions extract: ${sample}_sorted.bam.junc ..."
    regtools junctions extract \
    -s 0 -a 8 -m 50 -M 500000 \
    /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam \
    -o /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam.junc
    cat /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted.bam.junc | grep "^chr" | grep -v "^chrM" > /path/to/GEUV/BAM/${sample}_STAR_PB29/${sample}_sorted_filtered.bam.junc
done



# EGA and ImmVar and DICE and GEUV

touch /path/to/leafcutter/all_juncfiles.txt
ls /path/to/EGA/BAM/*_STAR_PB29/*_sorted_filtered.bam.junc >> /path/to/leafcutter/all_juncfiles.txt
ls /path/to/ImmVar/BAM/*_STAR_PB29/*_sorted_filtered.bam.junc >> /path/to/leafcutter/all_juncfiles.txt
ls /path/to/DICEBAM/*/*_STAR_PB29/*_sorted_filtered.bam.junc >> /path/to/leafcutter/all_juncfiles.txt
ls /path/to/GEUV/BAM/*_STAR_PB29/*_sorted_filtered.bam.junc >> /path/to/leafcutter/all_juncfiles.txt

nohup python tools/leafcutter/clustering/leafcutter_cluster_regtools.py \
-j /path/to/leafcutter/all_juncfiles.txt \
-r /path/to/leafcutter/ \
-o all \
> leafcutter_cluster_regtools.log &


.pyenv/versions/2.7.18/bin/python tools/leafcutter/scripts/prepare_phenotype_table_modified.py \
/path/to/leafcutter/all_perind.counts.gz -p 10

# make group files
## scripts/nanopore/make_groups_leafcutter.R

# ImmVar
for cell in CD4T_Activated MoDC_unstim MoDC_FLU MoDC_IFNb; do
    echo -e "ImmVar, ${cell}"
    ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
    --num_threads 30 \
    --output_prefix /path/to/leafcutter/leafcutter_ds_ImmVar_${cell} \
    --exon_file /path/to/tools/leafcutter/annotation_codes/PB29/PB29_all_exons.txt.gz \
    /path/to/leafcutter/all_perind_numers.counts.gz \
    /path/to/leafcutter/ImmVar_${cell}_groups.txt

    echo -e "ImmVar, ${cell}"
    ~/tools/leafcutter/leafviz/prepare_results_modified.R \
    /path/to/leafcutter/all_perind_numers.counts.gz \
    /path/to/leafcutter/leafcutter_ds_ImmVar_${cell}_cluster_significance.txt \
    /path/to/leafcutter/leafcutter_ds_ImmVar_${cell}_effect_sizes.txt \
    /path/to/tools/leafcutter/annotation_codes/PB29/PB29 \
    --meta_data_file /path/to/leafcutter/ImmVar_${cell}_groups.txt \
    --code ${cell}_vs_others \
    --num_threads 30 \
    --ntop_pca 500 \
    --output /path/to/leafcutter/leafcutter_ds_ImmVar_${cell}_leafviz.RData
done

# EGA
for cell in NS LPS Pam3CSK4 R848 IAV; do
    echo -e "EGA, ${cell}"
    ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
    --num_threads 30 \
    --output_prefix /path/to/leafcutter/leafcutter_ds_EGA_${cell} \
    --exon_file /path/to/tools/leafcutter/annotation_codes/PB29/PB29_all_exons.txt.gz \
    /path/to/leafcutter/all_perind_numers.counts.gz \
    /path/to/leafcutter/EGA_${cell}_groups.txt
    
    echo -e "EGA, ${cell}"
    ~/tools/leafcutter/leafviz/prepare_results_modified.R \
    /path/to/leafcutter/all_perind_numers.counts.gz \
    /path/to/leafcutter/leafcutter_ds_EGA_${cell}_cluster_significance.txt \
    /path/to/leafcutter/leafcutter_ds_EGA_${cell}_effect_sizes.txt \
    /path/to/tools/leafcutter/annotation_codes/PB29/PB29 \
    --meta_data_file /path/to/leafcutter/EGA_${cell}_groups.txt \
    --code ${cell}_vs_others \
    --num_threads 10 \
    --ntop_pca 500 \
    --output /path/to/leafcutter/leafcutter_ds_EGA_${cell}_leafviz.RData
done


# DICE
cat /path/to/DICE/cell_type.list | while read cell; do
    echo -e "DICE, ${cell}"
    ~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
    --num_threads 30 \
    --output_prefix /path/to/leafcutter/leafcutter_ds_DICE_${cell} \
    --exon_file /path/to/tools/leafcutter/annotation_codes/PB29/PB29_all_exons.txt.gz \
    /path/to/leafcutter/all_perind_numers.counts.gz \
    /path/to/leafcutter/DICE_${cell}_groups.txt
    
    echo -e "DICE, ${cell}"
    ~/tools/leafcutter/leafviz/prepare_results_modified.R \
    /path/to/leafcutter/all_perind_numers.counts.gz \
    /path/to/leafcutter/leafcutter_ds_DICE_${cell}_cluster_significance.txt \
    /path/to/leafcutter/leafcutter_ds_DICE_${cell}_effect_sizes.txt \
    /path/to/tools/leafcutter/annotation_codes/PB29/PB29 \
    --meta_data_file /path/to/leafcutter/DICE_${cell}_groups.txt \
    --code ${cell}_vs_others \
    --num_threads 30 \
    --ntop_pca 500 \
    --output /path/to/leafcutter/leafcutter_ds_DICE_${cell}_leafviz.RData
done

# GEUV
echo -e "GEUV"
~/tools/leafcutter/scripts/leafcutter_ds_modified.R \
--num_threads 30 \
--output_prefix /path/to/leafcutter/leafcutter_ds_GEUV \
--exon_file /path/to/tools/leafcutter/annotation_codes/PB29/PB29_all_exons.txt.gz \
/path/to/leafcutter/all_perind_numers.counts.gz \
/path/to/leafcutter/GEUV_CEU1YRI0_groups.txt

echo -e "ImmVar, ${cell}"
~/tools/leafcutter/leafviz/prepare_results_modified.R \
/path/to/leafcutter/all_perind_numers.counts.gz \
/path/to/leafcutter/leafcutter_ds_GEUV_cluster_significance.txt \
/path/to/leafcutter/leafcutter_ds_GEUV_effect_sizes.txt \
/path/to/tools/leafcutter/annotation_codes/PB29/PB29 \
--meta_data_file /path/to/leafcutter/GEUV_CEU1YRI0_groups.txt \
--code CEU_vs_YRI \
--num_threads 30 \
--ntop_pca 500 \
--output /path/to/leafcutter/leafcutter_ds_GEUV_leafviz.RData



