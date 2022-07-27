#!/bin/sh
#$ -S /bin/sh

# SUPPA

# Generate the ioe files:

## sr-RNAseq corrected
bio/python3/bin/python3 /path/to/tools/SUPPA-2.3/suppa.py generateEvents \
-i /path/to/GRCh38_gencode38_classification.filtered_lite.gtf \
-o /path/to/tools/SUPPA-2.3/SUPPA_supplementary_data/annotation/Homo_sapiens.GRCh38.PB29-all-sr-corrected.annotation.events \
-e SE SS MX RI FL -f ioe
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' /path/to/tools/SUPPA-2.3/SUPPA_supplementary_data/annotation/Homo_sapiens.GRCh38.PB29-all-sr-corrected.annotation.events*.ioe \
> /path/to/tools/SUPPA-2.3/SUPPA_supplementary_data/annotation/Homo_sapiens.GRCh38.PB29-all-sr-corrected.annotation.events.ioe
