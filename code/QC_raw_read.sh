#!/bin/sh
#$ -S /bin/sh

PBs=(
    "PB01"
    "PB02"
    "PB03"
    "PB04"
    "PB05"
    "PB06"
    "PB07"
    "PB08"
    "PB09"
    "PB10"
    "PB11"
    "PB12"
    "PB13"
    "PB14"
    "PB15"
    "PB16"
    "PB17"
    "PB18"
    "PB19"
    "PB20"
    "PB21_A"
    "PB21_B"
    "PB21_C"
    "PB22"
    "PB23"
    "PB24"
    "PB25"
    "PB26"
    "PB27"
    "PB28"
    "PB29"
)

cells=(
    "NaiveCD4"
	"Th1"
	"Th2"
	"Th17"
	"Tfh"
	"Fra1"
	"Fra2-aTreg"
	"Fra3"
	"LAG3Treg"
	"MemoryCD4"
	"Thx"
	"NaiveCD8"
	"CD8effector"
	"CD8centralmem"
	"CD8effectormem"
	"NaiveB"
	"unswmemoryB"
	"swmemoryB"
	"DNB"
	"plasmablast"
	"plasmacytoidDC_A"
	"plasmacytoidDC_B"
	"plasmacytoidDC_C"
	"myeloidDC"
	"NK"
	"monocyteCD16"
	"monocyteCD16minus"
	"nonclassicalMonocyte"
	"intermediateMonocyte"
	"PBMC"
	"Neutrophil"
)


for n in {0..30} ; do
echo "${PBs[n]}; ${cells[n]}"
mkdir --parents /path/to/${cells[n]}/NanoPlot

# raw reads
NanoPlot \
-t 10 \
--fastq /path/to/${PBs[n]}*/*/fastq_*.fastq.gz \
--title ${cells[n]}_fastq --N50 --tsv_stats -o /path/to/${cells[n]}/NanoPlot -p ${cells[n]}_fastq_NanoPlot \
--plots hex dot --legacy hex

# mapped reads
NanoPlot \
-t 20 \
--color red \
--bam /path/to/${cells[n]}/flair/GRCh38.aligned.bam \
--title ${cells[n]}_bam --N50 --tsv_stats -o /path/to/${cells[n]}/NanoPlot -p ${cells[n]}_bam_NanoPlot \
--plots hex dot --legacy hex
done

