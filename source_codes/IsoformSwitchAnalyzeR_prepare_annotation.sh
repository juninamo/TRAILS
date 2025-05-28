#!/bin/bash
#$ -S /bin/sh

# IsoformSwitchAnalyzeR

# prepare extract annotation
## Pfam
pfam_scan.pl -fasta /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_aa.fasta -dir tools/PfamScan/Pfam35.0 > /path/to/IsoformSwitchAnalyzeR/pfam_all.result
pfam_scan.pl -fasta /path/to/IsoformSwitchAnalyzeR/isoformSwitchAnalyzeR_isoform_AA_complete.fasta -dir tools/PfamScan/Pfam35.0 > /path/to/IsoformSwitchAnalyzeR/pfam.result

## CPAT
fread("/path/to/cpat/output.ORF_prob.best.tsv") %>%
  dplyr::mutate(`Coding Label` = ifelse(Coding_prob >= 0.364, "yes", "no"),
                `Sequence Name` = ifelse(!grepl("^ENST",seq_ID),tolower(seq_ID), seq_ID),
                `RNA size` = mRNA,
                `ORF size` = ORF,
                `Ficket Score` = Fickett,
                `Hexamer Score` = Hexamer,
                `Coding Probability` = Coding_prob) %>%
  dplyr::arrange(`Sequence Name`) %>%
  dplyr::mutate(`Data ID` = dplyr::row_number()) %>%
  as.data.frame() %>%
  .[,c("Data ID","Sequence Name","RNA size","ORF size","Ficket Score","Hexamer Score","Coding Probability","Coding Label")] %>%
  write.table(.,"/path/to/IsoformSwitchAnalyzeR/cpat_all.result",quote = FALSE,sep="\t",row.names = FALSE)

## IUPred
tool="SQANTI3"
path_to_aafasta="/path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_aa.fasta"
grep "^>" ${path_to_aafasta} | sed -e "s/^>//g" | awk '{print $1}' > read_id_${tool}.list
cat read_id_${tool}.list | while read read_id; do
grep -v "^#" /path/to/IUPred3_output/${tool}/${read_id}.txt | sed -e "1i# POS\tAMINO ACID\tIUPRED SCORE\tANCHOR SCORE" | sed -e "1i>${read_id}" | sed -e "1i# Nucleic Acids Research 2018, Submitted" | sed -e "1i# Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi" | sed -e "1i# IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding" | sed -e "1i################" > tmp/${read_id}_iupred.txt
done
find tmp -name "*_iupred.txt" | xargs cat > /path/to/IsoformSwitchAnalyzeR/iupred_all.result

## signalP
seqkit split -p 9 /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite_aa.fasta -O tmp/fasta
for id in {1..9}; do
echo $id
cd ~/tools/signalp-4.1
signalp \
-f summary \
~/tmp/fasta/GRCh38_gencode38_classification.filtered_lite_aa.part_00${id}.fasta \
> ~/tmp/fasta/GRCh38_gencode38_classification.filtered_lite_aa.part_00${id}.result
done
cd
rm ~/tmp/sigp_*.txt
rm ~/doc*
for id in {1..9}; do
awk '/^# Measure/{f="doc."++d} f{print > f}' ~/tmp/fasta/GRCh38_gencode38_classification.filtered_lite_aa.part_00${id}.result
ls doc.* | while read line; do
seq_id=$(grep -e "^Name=" $line | awk '{print $1}' | sed -e "s/Name=//g")
cat $line | sed -e "1i>${seq_id}" | grep -v "# SignalP-4.1 euk predictions" > ~/tmp/sigp_${seq_id}.txt
done
rm ~/doc*
done
find tmp -name "sigp_*.txt" | xargs cat | sed -e "1i# SignalP-4.1 euk predictions" > /path/to/IsoformSwitchAnalyzeR/signalp_all.result
rm ~/tmp/sigp_*.txt
