#!/bin/sh
#$ -S /bin/sh

# SQUANTI3

# CAGE
cat \
/path/to/reference/refTSS/v3.3/human/refTSS_v3.3_human_coordinate.hg38.bed \
/path/to/reference/FANTOM5/TSS_classifier/TSS_human_hg38_relaxed_strict_TSS.bed \
| LANG=C sort -V -k 1,1 -k 2,2n > /path/to/tools/SQANTI3-4.2/data/ref_TSS_annotation/merged_refTSS-v3.3_TSS-classifier-relaxed-strict_hg38.bed


python /path/to/tools/SQANTI3-4.2/sqanti3_qc.py \
/path/to/PB29_common/GRCh38_all-sr-correct.isoforms.gtf \
/path/to/reference/GENCODE/GRCh38/release_38/gencode.v38.annotation.gtf \
/path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
--min_ref_len 200 \
--aligner_choice minimap2 \
--coverage /path/to/tools/SQANTI3-4.2/data/sr_RNAseq/all-sr_sum_SJ_thresh3.out.tab \
--sites ATAC,GCAG,GTAG \
--window 20 \
--genename \
--cage_peak /path/to/tools/SQANTI3-4.2/data/ref_TSS_annotation/merged_refTSS-v3.3_TSS-classifier-relaxed-strict_hg38.bed \
--polyA_motif_list /path/to/tools/SQANTI3-4.2/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
--polyA_peak /path/to/tools/SQANTI3-4.2/data/polyA_site/atlas.clusters.2.0.GRCh38.96.bed \
--output GRCh38_gencode38 \
--dir /path/to/PB29_common/SQANTI3_output \
--cpus 4 \
--saturation \
--report both \
--isoAnnotLite \
--gff3 /path/to/tools/tappAS-1.0.7/resources/gffs/Homo_sapiens_GRCh38_Ensembl_86.gff3

python /path/to/tools/SQANTI3-4.2/sqanti3_RulesFilter.py \
/path/to/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.txt \
/path/to/PB29_common/SQANTI3_output/GRCh38_gencode38_corrected.fasta \
/path/to/PB29_common/SQANTI3_output/GRCh38_gencode38_corrected.gtf \
--intrapriming 0.6 \
--max_dist_to_known_end 50 \
--min_cov 3 \
--saturation \
--report both

# generate transcript sequence
gffread \
/path/to/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
-g /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-w /path/to/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fa


