#!/bin/sh
#$ -S /bin/sh

# Prjibelski, A.D., Mikheenko, A., Joglekar, A. et al. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-022-01565-y

# we used "--report_novel_unspliced true" to compare to FLAIR pipeline, which also reports mono-exon
/path/to/IsoQuant-3.2.0/isoquant.py \
--data_typ nanopore \
--bam_list /path/to/bam_list.txt \
--reference /path/to/reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
--stranded none \
--model_construction_strategy sensitive_ont \
--report_novel_unspliced true \
--genedb /path/to/gencode.v38.annotation.gtf \
--complete_genedb \
--sqanti_output \
--threads 16 \
-o /path/to/isoquant