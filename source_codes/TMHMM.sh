#!/bin/bash
##$ -S /bin/sh

tools/TMHMM/tmhmm-2.0c/bin/tmhmm -short /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_aa.fasta > /path/to/TMHMM_output/SQANTI3_output.txt
tools/TMHMM/tmhmm-2.0c/bin/tmhmm /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_aa.fasta > /path/to/TMHMM_output/SQANTI3_output_long.txt
