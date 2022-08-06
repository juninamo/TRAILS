#!/bin/bash
#$ -S /bin/sh

GEO_id="GSE61742"
sample_list=$(cat /path/to/riboseq/${GEO_id}/SRR_Acc_List.txt | tr "\n" " ")

#============================================
# step 0: download
#============================================

for ribo_id in $sample_list; do
mkdir -p /path/to/riboseq/${GEO_id}/dump/${ribo_id}
if [ -e "/path/to/riboseq/${GEO_id}/dump/${ribo_id}/${ribo_id}.fastq.gz" ]; then
echo -e "${ribo_id}.fastq.gz exists"
else
echo -e "${ribo_id}.fastq.gz not exists"
DIR="/path/to/riboseq/${GEO_id}/dump/${ribo_id}"
tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump-orig.2.11.2 -O ${DIR} --split-3 -e 4 -p ${ribo_id}
bin/pigz -p 4 ${DIR}/${ribo_id}.fastq
fi
done

#============================================
# step 1: trim adaptor sequence
#============================================

trim_out="/path/to/riboseq/trim_galore/${GEO_id}/"
sample_list=$(cat /path/to/riboseq/${GEO_id}/SRR_Acc_List.txt | tr "\n" " ")

for ribo_id in $sample_list; do
riboseq_fq="/path/to/riboseq/${GEO_id}/dump/${ribo_id}/${ribo_id}.fastq.gz"
if [ -e "${trim_out}/${ribo_id}_trimmed.fq.gz" ]; then
echo -e "${ribo_id}_trimmed.fq.gz exists"
else
trim_galore -q 20 -o ${trim_out} ${riboseq_fq}
fi
done

#============================================
# step 2: filter rrna
#============================================

# Ribo-seq: to remove rRNA
## https://imamachi-n.hatenablog.com/entry/2017/03/31/215218
### Human 5S DNA
#### https://www.ncbi.nlm.nih.gov/nuccore/23898
### Human ribosomal DNA complete repeating unit
#### https://www.ncbi.nlm.nih.gov/nuccore/555853
mkdir -p STAR/rRNA_contam
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir STAR/rRNA_contam \
--genomeFastaFiles reference/fasta/contam_Ribosomal_RNA.fa

nproc=8 # threads
rrna_idx="STAR/rRNA_contam"
trim_out="/path/to/riboseq/trim_galore/${GEO_id}/"
star_ribo_out="/path/to/riboseq/STAR/${GEO_id}/"
sample_list=$(cat /path/to/riboseq/${GEO_id}/SRR_Acc_List.txt | tr "\n" " ")

for ribo_id in $sample_list; do
if [ -e "${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastq" ]; then
echo -e "${ribo_id}_rm_repbase_rrna.fastq exists"
else
STAR \
--runMode alignReads \
--runThreadN $nproc \
--genomeDir $rrna_idx \
--readFilesIn ${trim_out}${ribo_id}_trimmed.fq.gz \
--readFilesCommand gunzip -c \
--outSAMunmapped Within \
--outFilterMultimapNmax 30 \
--outFilterMultimapScoreRange 1 \
--outFileNamePrefix ${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastq \
--outSAMattributes All \
--outStd BAM_Unsorted \
--outSAMtype BAM Unsorted \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10 \
--alignEndsType EndToEnd > ${star_ribo_out}${ribo_id}_repbase_rrna_comtam.bam
mv ${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastqUnmapped.out.mate1 ${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastq
fi
done

#============================================
# step 3: mapping to genome
#============================================

# Ribo-seq: mapping to genome
mkdir --parents STAR/genome_GRCh38_PB29
STAR --runMode genomeGenerate \
--genomeDir STAR/genome_GRCh38_PB29/ \
--genomeFastaFiles reference/GENCODE/GRCh38/release_38/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.gtf \
--sjdbOverhang 100 \
--runThreadN 35

nproc=8 # threads
# ribo_id="SRR1585486"
star_ribo_out="/path/to/riboseq/STAR/${GEO_id}/"
sample_list=$(cat /path/to/riboseq/${GEO_id}/SRR_Acc_List.txt | tr "\n" " ")

for ribo_id in $sample_list; do
if [ -e "${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_sorted_PB29.bam" ]; then
echo -e "${ribo_id}_rm_repbase_rrna_genome_sorted_PB29.bam exists"
else
STAR \
--runMode alignReads \
--runThreadN $nproc \
--genomeDir STAR/genome_GRCh38_PB29 \
--readFilesIn ${star_ribo_out}${ribo_id}_rm_repbase_rrna.fastq \
--outFilterMultimapNmax 8 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
--outFilterMismatchNmax 4 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFileNamePrefix ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_ \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 32000000000 \
--outFilterType BySJout \
--outReadsUnmapped Fastx

/path/to/samtools sort -@ 8 ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_Aligned.sortedByCoord.out.bam -o ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_sorted_PB29.bam
/path/to/samtools index -@ 8 ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_sorted_PB29.bam
rm ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_Aligned.sortedByCoord.out.bam
fi
done

#============================================
# step 4: quality check of riboseq and identify TIS by ribotish
#============================================

mkdir --parents /path/to/riboseq/ribotish/${GEO_id}
star_ribo_out="/path/to/riboseq/STAR/${GEO_id}/"
sample_list=$(cat /path/to/riboseq/${GEO_id}/SRR_Acc_List.txt | tr "\n" " ")

for ribo_id in $sample_list; do
if [ -e "/path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.pdf" ]; then
echo -e "/path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.pdf exists"
else
echo -e "ribotish quality: ${ribo_id} ..."
ribotish quality -p 30 -b ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_sorted_PB29.bam -g /data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.utr.start.stop.sorted.gtf --geneformat gtf -o /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.txt -f /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.pdf -r /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.para.py --verbose > /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.log 2>&1
fi
if [ -e "/path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_pred.txt" ]; then
echo -e "/path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_pred.txt exists"
else
echo -e "ribotish predict: ${ribo_id} ..."
ribotish predict -p 30 -b ${star_ribo_out}${ribo_id}_rm_repbase_rrna_genome_sorted_PB29.bam -g /data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.utr.start.stop.sorted.gtf --geneformat gtf -f reference/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --ribopara /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_qual.para.py --seq --aaseq --blocks --longest --alt --altcodons AAA,AAG,AAT,AAC,ACG,ACT,ACC,ACA,AGG,AGT,AGC,AGA,ATC,ATA,ATT,CTG,CAT,TAT,TAC,TTG,TGG,TGT,TGC,GTG,GCT,GCC,GCA,GCG,GGT,GGC,GGA,GGG,GAT,GAC,GTT,GTA,GTC,GTG -o /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_pred.txt --verbose > /path/to/riboseq/ribotish/${GEO_id}/${ribo_id}_pred.log 2>&1
fi
done

# OUTPUT
## The output is a txt file all possible ORF results that fit the thresholds. Some of the columns are:
### GenomePos:	Genome position and strand of TIS site, 0 based, half open
### Start:	TIS of the ORF on transcript
### Stop:	3' end of stop codon on transcript
### TisType:	Relative position of this TIS to annotated ORF of the transcript. 'Novel' if no ORF annotation. ':Known' means the TIS is annotated in another transcript. ':CDSOverlap' means the ORF overlaps with annotated CDS in another transcript in the same reading frame.
### TISGroup:	Group of the transcript for TIS background estimation
### TISCount:	Number of reads with P-site at TIS site
### TISPvalue:	One tailed negative binomial test p-value for TISCount (TIS test)
### RiboPvalue:	One tailed rank sum test p-value for regular riboseq frame bias inside ORF (frame test)
### RiboPStatus:	For all ORFs sharing same stop codon, 'T' means top (best) p-value, 'L' means local best p-value, 'N' means other. All 'N' in `-i` or `--longest` mode.
### FisherPvalue:	Combination of TIS and Ribo p-values using Fisher's method
### TISQvalue:	BH correction q-value of TIS test
### RiboQvalue:	BH correction q-value of frame test
### FisherQvalue:	BH correction q-value of Fisher's p-value
### AALen:	Amino acid length of the ORF

#============================================
# step 5: identify uORFs by modified scripts of uORF-tools, count reads amd calculate uORF/CDS ratio of individual sample
#============================================

mkdir --parents /path/to/riboseq/uORF/${GEO_id}
python /path/to/uORF-Tools/scripts/longest_orf_transcript_modified.py -a /data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.utr.start.stop.sorted.gtf -o /path/to/riboseq/uORF/longest_protein_coding_transcripts.gtf

cat /path/to/riboseq/ribotish/${GEO_id}/SRR*_pred.txt | head -n 1 > hdr.txt
cat /path/to/riboseq/ribotish/${GEO_id}/SRR*_pred.txt | grep -v "^Gid" > all.txt
cat hdr.txt all.txt > /path/to/riboseq/ribotish/${GEO_id}/merged_pred.txt

echo -e "riboMerge ..."
bed="/path/to/riboseq/uORF/${GEO_id}/merged_uORFs.bed"
csv="/path/to/riboseq/uORF/${GEO_id}/merged_uORFs.csv"
python /path/to/uORF-Tools/scripts/ribo_merge.py /path/to/riboseq/ribotish/${GEO_id}/merged_pred.txt --min_length 1 --max_length 400 --output_csv_filepath ${csv} --output_bed_filepath ${bed}

echo -e "size_factors ..."
star_ribo_out="/path/to/riboseq/STAR/${GEO_id}/"
Rscript /path/to/uORF-Tools/scripts/generate_size_factors_modified.R -b ${star_ribo_out} -a /data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.utr.start.stop.sorted.gtf -t /path/to/riboseq/uORF/${GEO_id}/yoruba_sample.tsv -s /path/to/riboseq/uORF/${GEO_id}/sfactors_lprot.csv

echo -e "count reads ..."
norm="/path/to/riboseq/uORF/${GEO_id}/ribo_norm_CDS_reads.csv"
raw="/path/to/riboseq/uORF/${GEO_id}/ribo_raw_CDS_reads.csv"
Rscript /path/to/uORF-Tools/scripts/generate_ribo_counts_CDS_modified.R -b ${star_ribo_out} -a /data03/inamo/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.cds.utr.start.stop.sorted.gtf -s /path/to/riboseq/uORF/${GEO_id}/sfactors_lprot.csv -t /path/to/riboseq/uORF/${GEO_id}/yoruba_sample.tsv -n ${norm} -r ${raw}
norm="/path/to/riboseq/uORF/${GEO_id}/ribo_norm_uORFs_reads.csv"
raw="/path/to/riboseq/uORF/${GEO_id}/ribo_raw_uORFs_reads.csv"
Rscript /path/to/uORF-Tools/scripts/generate_ribo_counts_uORFs_modified.R -b ${star_ribo_out} -a /path/to/riboseq/uORF/${GEO_id}/merged_uORFs.bed -s /path/to/riboseq/uORF/${GEO_id}/sfactors_lprot.csv -t /path/to/riboseq/uORF/${GEO_id}/yoruba_sample.tsv -n ${norm} -r ${raw}

echo -e "ribo_changes ..."
uorf="/path/to/riboseq/uORF/${GEO_id}/ribo_norm_uORFs_reads.csv"
orf="/path/to/riboseq/uORF/${GEO_id}/ribo_norm_CDS_reads.csv"
frac="/path/to/riboseq/uORF/${GEO_id}/ribo_change.csv"
python /path/to/uORF-Tools/scripts/ribo_changes.py --uORF_reads ${uorf} --ORF_read ${orf} --changes_output ${frac}

echo -e "final table ..."
annotation="/path/to/riboseq/uORF/${GEO_id}/merged_uORFs.csv"
uORFreads="/path/to/riboseq/uORF/${GEO_id}/ribo_norm_uORFs_reads.csv"
cdsreads="/path/to/riboseq/uORF/${GEO_id}/ribo_norm_CDS_reads.csv"
output="/path/to/riboseq/uORF/${GEO_id}/uORFs_regulation.tsv"
python /path/to/uORF-Tools/scripts/final_table.py --uORF_reads ${uORFreads} --ORF_reads ${cdsreads} --uORF_annotation ${annotation} --output_csv_filepath ${output}

