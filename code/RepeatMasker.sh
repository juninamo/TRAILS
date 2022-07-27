#!/bin/sh
#$ -S /bin/sh

# RepeatMasker
## https://rpubs.com/nishikosh/RepeatMasker_output

RepeatMasker -e hmmer -pa 20 /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fasta -species human -s -nolow -gff -html -ace -source -dir /path/to/rmout \
> /path/to/`date "+%Y%m%d"`_RepeatMasker.log 2>&1

perl /path/to/tools/RepeatMasker/util/rmOutToGFF3.pl /path/to/rmout/GRCh38_gencode38_classification.filtered_lite.fasta.out > /path/to/rmout/GRCh38_gencode38_classification.filtered_lite.fasta.out.gff3
perl /path/to/tools/RepeatMasker/util/rmOut2Fasta.pl -fasta /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fasta -out /path/to/rmout/GRCh38_gencode38_classification.filtered_lite.fasta.out > /path/to/rmout/GRCh38_gencode38_classification.filtered_lite.fa
