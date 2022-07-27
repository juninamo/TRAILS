#!/bin/bash
#$ -S /bin/sh

# Nature 505, 706–709 (2014).
## Methods Summary
### Sample preparation and structure probing for human renatured RNAs;
#### Human lymphoblastoid cell lines GM12878, GM12891 and GM12892 were obtained from Coriell. Total RNA was isolated using TRIzol reagent (Invitrogen) and polyA selected as described previously2. Two micrograms of Poly(A)+ RNA was structure probed at 37 °C using RNase V1 (Life Technologies, final concentration of 10−5 units per μl) or S1 nuclease (Fermentas, final concentration of 0.4 units per μl) at 37 °C for 15 min.
### Sample preparation and structure probing for human native deproteinized RNAs;
#### GM12878 cells were lysed in lysis buffer (150 mM NaCl, 10 mM MgCl2, 1% NP40, 0.1% SDS, 0.25% Na deoxycholate, Tris pH 7.4) on ice for 30 min. The lysate was deproteinized by phenol chloroform extractions. Total RNA (1 μg per 90 μl) was incubated in 1 × RNA structure buffer at 37 °C for 15 min and structure probed using RNase V1 (final concentration of 2 × 10−5 units per μl) and S1 nuclease (final concentration of 0.2 units per μl) at 37 °C for 15 min.
### Library construction and analysis;
#### The structure probed RNA was cloned using Ambion RNA-Seq Library Construction Kit (Life Technologies)2, and sequenced using Illumina Hi-seq. The reads were trimmed and mapped to UCSC RefSeq and the Gencode v12 databases (hg19 assembly) using the software Bowtie2 (ref. 17). Double (V1) and single-stranded reads (S1) for each sequencing sample were normalized by sequencing depth.
### RiboSNitch analysis;
#### Data normalization for each sample was performed by calculating standard deviation (s.d.) for each transcript and dividing the PARS score per base by the s.d. of that transcript. We defined a structure difference of the ith base of transcript j between conditions m and n in this formula, where PARS represents the normalized PARS score, abs represents absolute value, and k represents the kth base of the transcript:
## Online Methods
### Sample preparation for renatured RNA structure probing;
#### Human lymphoblastoid cell lines GM12878, GM12891 and GM12892 were obtained from Coriell. Total RNA was isolated from lymphoblastoid cells using TRIzol reagent (Invitrogen). Poly(A)+ RNA was obtained by purifying twice using the MicroPoly(A)Purist kit (Life Technologies). The Tetrahymena ribozyme RNA was in vitro transcribed using the T7 RiboMax Large-scale RNA production system (Promega) and added into 2 μg of poly(A)+ RNA (1% by mole) for structure probing and library construction.
### Structure probing of renatured poly(A)+ RNA;
#### Two micrograms of Poly(A)+ RNA in 160 μl of nuclease free water was heated at 90 °C for 2 min and snap-cooled on ice for 2 min. Twenty microlitres of 10 × RNA structure buffer (150 mM NaCl, 10 mM MgCl2, Tris, pH 7.4) was added to the RNA and the RNA was slowly warmed up to 37 °C over 20 min. The RNA was then incubated at 37 °C for 15 min and structure probed independently using RNase V1 (Life Technologies, final concentration of 10−5 units per μl) or S1 nuclease (Fermentas, final concentration of 0.4 units per μl) at 37 °C for 15 min. The cleavage reactions were inactivated using phenol chloroform extraction.
### Structure probing and ribosomal RNA depletion for native deproteinized RNA structure probing;
#### GM12878 cells were lysed in lysis buffer (150 mM NaCl, 10 mM MgCl2, 1% NP40, 0.1% SDS, 0.25% Na deoxycholate, Tris, pH 7.4) on ice for 30 min. The chromatin pellet was removed by centrifugation at 16,000g for 10 min at 4 °C. The lysate was deproteinized by passing through two phenol followed by one chloroform extractions. The concentration of RNA in the deproteinized lysate was measured using the Qubit fluorometer (Invitrogen). We diluted the RNA to a concentration of 1 μg per 90 μl using 1 × RNA structure buffer (150 mM NaCl, 10 mM MgCl2, Tris, pH 7.4) and incubated the RNA at 37 °C for 15 min. The native deproteinized RNA was structure probed independently using RNase V1 (final concentration of 2 × 10−5 units per μl) and S1 nuclease (final concentration of 0.2 units per μl) at 37 °C for 15 min.
#### To compare structural differences between renatured and native deproteinized RNAs, we independently prepared an RNA sample that was similarly lysed and deproteinized. After removal of proteins, we ethanol precipitated the RNA and dissolved it in nuclease free water. We diluted the RNA to a concentration of 1 μg per 80 μl in water and heated the RNA at 90 °C for 2 min before snap-cooling the RNA on ice. We added 10 × RNA structure buffer and renatured the RNA by incubating it at 37 °C for 15 min and performed structure probing similarly as in native deproteinized RNAs.
#### The cleavage reactions were inactivated using phenol chloroform extraction and DNase treated before undergoing ribosomal RNA depletion using Ribo-Zero Ribosomal RNA removal kit (Epicentre).
### Validation of riboSNitches by manual footprinting;
#### We cloned approximately 200 nucleotide fragments of both alleles of MRPS21, WSB1, HLA-DRB1, HLA-DQA1, hnRNP-AB, HLA-DRA, LDHA, XRCC5 and FNBP1 from GM12878, GM12891 and GM12892 using a forward-T7-gene-specific primer and a reverse-gene-specific primer. All constructs were confirmed by sequencing using capillary electrophoresis. DNA from each of the different clones was then in vitro transcribed into RNA using MegaScript Kit from Ambion, following manufacturer’s instructions.
#### Two picomoles of each RNA is heated at 90 °C for 2 min and chilled on ice for 2 min. 3.33 × RNA folding mix (333 mM HEPES, pH 8.0, 20 mM MgCl2, 333 mM NaCl) was then added to the RNA and the RNA was allowed to fold slowly to 37 °C over 20 min. The RNA was then structure probed with either DMS (final concentration of 100 mM) or 2-methylnicotinic acid imidazolide (NAI) (final concentration of 100 mM)16 at 37 °C for 20 min or structure probed with S1 nuclease (final concentration of 0.4 units per µl) or RNase V1 (final concentration of 0.0001 units per µl) at 37 °C for 15 min. The DMS structure probed samples were quenched using 2-mercaptoethanol before phenol chloroform extraction. The NAI and nuclease treated samples were phenol chloroform extracted directly after structure probing. The structure probed RNA was then recovered through ethanol precipitation. The RNA structure modification/cleavage sites were then read out using a radiolabelled RT primer by running onto denaturing PAGE gel as described previously18.
### Library construction;
#### The structure-probed RNA was fragmented at 95 °C using alkaline hydrolysis buffer (50 mM Sodium Carbonate, pH 9.2, 1 mM EDTA) for 3.5 min. The fragmented RNA was then ligated to 5′ and 3′ adapters in the Ambion RNA-Seq Library Construction Kit (Life Technologies). The RNA was then treated with Antarctic phosphatase (NEB) to remove 3′ phosphates before re-ligating using adapters in the Ambion RNA-Seq Library Construction Kit (Life Technologies). The RNA was reverse-transcribed using 4 μl of the RT primer provided in the Ambion RNA-Seq Library Construction Kit and polymerase chain reaction (PCR)-amplified following the manufacturer’s instructions. We performed 18 cycles of PCR to generate the complementary DNA library.
### Illumina sequencing and mapping;
#### We performed paired end sequencing on Illumina’s Hi-Seq sequencer and obtained approximately 400-million reads for each paired end lane in an RNase V1 or S1 nuclease library. Obtained raw reads were truncated to 50 bases, (51 bases from the 3′ end were trimmed). Trimmed reads were mapped to the human transcriptome, which consists of non-redundant transcripts from UCSC RefSeq and the Gencode v12 databases (hg19 assembly), using the software Bowtie2 (ref. 17). We allowed up to one mismatch per seed during alignment, and only included reads with perfect mapping or with Bowtie2 reported mismatches on positions annotated as SNVs in genetically modified cells. We obtained 166- to 212-million mapped reads for an RNase V1 or S1 nuclease sample.
### PARS-score calculation;
#### After the raw reads were mapped to the transcriptome, we calculated the number of double-stranded reads and single-stranded reads that initiated on each base on an RNA. The number of double (V1) and single stranded reads (S1) for each sequencing sample were then normalized by sequencing depth. For a transcript with N bases in total, the PARS score of its ith base was defined by the following formula where V1 and S1 are normalized V1 and S1 scores, respectively. A small number 5 was added to reduce the potential over-estimating of structural signals of bases with low coverage:
#### To identify structural changes caused by SNVs, we applied a 5-base average on the normalized V1 and S1 scores to smoothing the nearby bases’ structural signals; therefore, the PARS score is defined as:
#### Bases with both high V1 and S1 scores, and transcripts with multiple conformations
#### Bases with both strong single- and double-strand signals are potentially present in multiple conformations. We first normalized all bases with detectable S1 or V1 counts by their sequencing depth. We then calculated an S1 ratio and a V1 ratio by normalizing S1 (and V1) counts to the transcript abundance. S1 and V1 ratios indicate the relative strength of single and double signals respectively. We then ranked all the bases by their S1 ratio and V1 ratio independently, and used the top one-million S1 ratio bases and the top one-million V1 ratio bases as high S1 ratio bases and high V1 ratio bases, respectively. We defined a base as being in multiple conformations if the base has both high S1 and high V1 ratios. If a transcript contains more than five multi-confirmation bases, this transcript is defined as a multi-confirmation transcript.
### V1 replicates correlation analysis;
#### Pearson correlation of RNase V1 replicates on GM12878 was performed using a parsV1 score (a value that uses the V1 score only to represent secondary structure) defined as:

GEO_id="GSE50676"
echo $GEO_id

for rna_id in SRR972706 SRR972707 SRR972708 SRR972709 SRR972710 SRR972711 SRR972712 SRR972713 SRR972714 SRR972715 SRR972716 SRR972717; do
mkdir -p /path/to/PARS/${GEO_id}/dump/${rna_id}
DIR="/path/to/PARS/${GEO_id}/dump/${rna_id}"
tools/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump-orig.2.11.2 -O ${DIR} --split-3 -e 4 -p ${rna_id}
bin/pigz -p 4 ${DIR}/${rna_id}.fastq
bin/pigz -p 4 ${DIR}/${rna_id}_1.fastq
bin/pigz -p 4 ${DIR}/${rna_id}_2.fastq
done

for rna_id in SRR972706 SRR972707 SRR972708 SRR972709 SRR972710 SRR972711 SRR972712 SRR972713 SRR972714 SRR972715 SRR972716 SRR972717; do
DIR="/path/to/PARS/${GEO_id}/dump/${rna_id}"
trim_galore -q 4 -o /path/to/PARS/${GEO_id}/trim_galore/ --paired /path/to/PARS/${GEO_id}/dump/${rna_id}/${rna_id}_1.fastq.gz /path/to/PARS/${GEO_id}/dump/${rna_id}/${rna_id}_2.fastq.gz
done


# rRNA depleted RNA and rRNA depleted RNA, native deproteinized replicate 1+2
# → mapping to isoforms which express in all LCL samples (GEUV, StringTie2+kallisto), nonredundunt (most expressed isoforms in associated gene) and predicted-coding by SQANTI3
## --outFilterMultimapNmax ; max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
## --alignSJoverhangMin ; minimum overhang for unannotated junctions
## --alignSJDBoverhangMin ; minimum overhang for annotated junctions
## --sjdbScore ;  extra alignment score for alignmets that cross database junctions
## --outFilterMismatchNmax ; maximum number of mismatches per pair, large number switches off this filter
## --alignIntronMin ; minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion
## --alignIntronMax ; maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins
## --alignMatesGapMax ; maximum gap between two mates, if 0, max intron gap will be determined by (2ˆwinBinNbits)*winAnchorDistNbins
### Interpreting the Log.final.out file - % of reads unmapped: too short? [https://groups.google.com/g/rna-star/c/jmszy_UjyKY/m/2drCEkG3CQAJ]
## --outFilterScoreMinOverLread ; real: same as outFilterScoreMin[int: alignment will be output only if its score is higher than or equal to this value.], but normalized to read length (sum of mates’lengths for paired-end reads)
## --outFilterMatchNminOverLread ; real: sam as outFilterMatchNmin[int: alignment will be output only if the number of matched bases is higher than or equal to this value], but normalized to the read length (sum of mates’ lengths for paired-end reads).
for rna_id in SRR972712 SRR972713 SRR972714 SRR972715 SRR972716 SRR972717; do
mkdir -p /path/to/PARS/${GEO_id}/STAR/${rna_id}_STAR_PB29_GEUVexp_nonredundunt
/path/to/STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir tools/STAR/genome_GRCh38_PB29_transcript_GEUVexp_nonredundunt \
--readFilesIn /path/to/PARS/${GEO_id}/trim_galore/${rna_id}_1_val_1.fq.gz /path/to/PARS/${GEO_id}/trim_galore/${rna_id}_2_val_2.fq.gz \
--readFilesCommand gunzip -c \
--outFilterMultimapNmax 8 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
--outFilterMismatchNmax 4 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFileNamePrefix /path/to/PARS/${GEO_id}/STAR/${rna_id}_STAR_PB29_GEUVexp_nonredundunt/ \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 32000000000 \
--outFilterType BySJout \
--outReadsUnmapped Fastx
/path/to/samtools sort -@ 8 /path/to/PARS/${GEO_id}/STAR/${rna_id}_STAR_PB29_GEUVexp_nonredundunt/Aligned.sortedByCoord.out.bam -o /path/to/PARS/${GEO_id}/STAR/${rna_id}_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam
/path/to/samtools index -@ 8 /path/to/PARS/${GEO_id}/STAR/${rna_id}_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam
rm /path/to/PARS/${GEO_id}/STAR/${rna_id}_STAR_PB29_GEUVexp_nonredundunt/Aligned.sortedByCoord.out.bam
done
rm /path/to/PARS/${GEO_id}/STAR/*_STAR_PB29_GEUVexp_nonredundunt/Unmapped.out.mate1
rm /path/to/PARS/${GEO_id}/STAR/*_STAR_PB29_GEUVexp_nonredundunt/Unmapped.out.mate2
mv /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam
mv /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam
mv /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam
mv /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam
mv /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam
mv /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam
mv /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam.bai
mv /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam.bai
mv /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam.bai
mv /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam.bai
mv /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam.bai
mv /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam.bai



# Number of reads mapped to each position in each transcript
## -a  output all positions (including zero depth)
## -a -a (or -aa)  output absolutely all positions, including unused ref. sequences
## -l <int> read length threshold (ignore reads shorter than <int>) [0]
## -d/-m <int> maximum coverage depth [8000]. If 0, depth is set to the maximum integer value, effectively removing any depth limit.
## -q <int> base quality threshold [0]
## -Q <int> mapping quality threshold [0]
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam > data03/inamo/PARS/GSE50676/STAR/SRR972717_STAR_PB29/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.depth.tab
/path/to/samtools depth -q 20 -Q 255 -aa -d 0 /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.depth.tab
# Number of reads mapped to each transcript
## reference sequence name, sequence length, mapped read-segments and unmapped read-segments
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.idxstats.tab
/path/to/samtools idxstats -@ 8 /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam > /path/to/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.idxstats.tab
mv /data03/inamo/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /data03/inamo/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam
mv /data03/inamo/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /data03/inamo/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam
mv /data03/inamo/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /data03/inamo/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam
mv /data03/inamo/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /data03/inamo/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam
mv /data03/inamo/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /data03/inamo/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam
mv /data03/inamo/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam /data03/inamo/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam
mv /data03/inamo/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /data03/inamo/PARS/GSE50676/STAR/SRR972712_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.bam.bai
mv /data03/inamo/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /data03/inamo/PARS/GSE50676/STAR/SRR972713_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.bam.bai
mv /data03/inamo/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /data03/inamo/PARS/GSE50676/STAR/SRR972714_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.bam.bai
mv /data03/inamo/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /data03/inamo/PARS/GSE50676/STAR/SRR972715_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.bam.bai
mv /data03/inamo/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /data03/inamo/PARS/GSE50676/STAR/SRR972716_STAR_PB29_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.bam.bai
mv /data03/inamo/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/transcriptome_sorted.bam.bai /data03/inamo/PARS/GSE50676/STAR/SRR972717_STAR_PB29_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.bam.bai




# https://github.com/dincarnato/RNAFramework
# https://rnaframework-docs.readthedocs.io/en/latest/

# rf-count
tools/RNAFramework/rf-count \
-p 30 -wt 30 \
-s /path/to/samtools -r \
-o /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt -ow \
-f /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.filtered.fasta \
/path/to/PARS/GSE50676/STAR/SRR9727*_STAR_PB29_GEUVexp_nonredundunt/GM12878_*_rRNAdepANDrenatured_transcriptome_sorted.bam

tools/RNAFramework/rf-count \
-p 30 -wt 30 \
-s /path/to/samtools -r \
-o /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt -ow \
-f /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.filtered.fasta \
/path/to/PARS/GSE50676/STAR/SRR9727*_STAR_PB29_GEUVexp_nonredundunt/GM12878_*_rRNAdep_rep1_transcriptome_sorted.bam

tools/RNAFramework/rf-count \
-p 30 -wt 30 \
-s /path/to/samtools -r \
-o /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt -ow \
-f /path/to/nanopore/RNA/210931_PB29_fastq/PB29_common/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.filtered.fasta \
/path/to/PARS/GSE50676/STAR/SRR9727*_STAR_PB29_GEUVexp_nonredundunt/GM12878_*_rRNAdep_rep2_transcriptome_sorted.bam


# rf-rctools
## column 1: Transcript ID
## column 2: Transcript sequence
## column 3: Number of per-base RT-stops (or mutations)
## column 4: Per-base coverage
tools/RNAFramework/rf-rctools view /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.tab
tools/RNAFramework/rf-rctools view /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.tab
tools/RNAFramework/rf-rctools stats /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.rc > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdepANDrenatured_transcriptome_sorted.stats
tools/RNAFramework/rf-rctools stats /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.rc > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdepANDrenatured_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdepANDrenatured_transcriptome_sorted.stats
tools/RNAFramework/rf-rctools view /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.tab
tools/RNAFramework/rf-rctools view /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.tab
tools/RNAFramework/rf-rctools stats /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.rc > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep1_transcriptome_sorted.stats
tools/RNAFramework/rf-rctools stats /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.rc > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep1_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep1_transcriptome_sorted.stats
tools/RNAFramework/rf-rctools view /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.tab
tools/RNAFramework/rf-rctools view /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.rc | tr "\n\n" "\t\t"  | sed -e "s/\t\t/\n/g" > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.tab
tools/RNAFramework/rf-rctools stats /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.rc > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_V1_rRNAdep_rep2_transcriptome_sorted.stats
tools/RNAFramework/rf-rctools stats /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.rc > /path/to/PARS/GSE50676/RNA_Framework/GM12878_rRNAdep_rep2_count_GEUVexp_nonredundunt/GM12878_S1_rRNAdep_rep2_transcriptome_sorted.stats



