#!/bin/bash
##$ -S /bin/sh

cd tools/signalp-5.0b/bin
signalp \
-format short \
-gff3 \
-fasta /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_aa.fasta \
-prefix /path/to/singalP_output/SQANTI3_signalp_output
