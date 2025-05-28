#!/bin/bash
#$ -S /bin/sh

# entire sequence
mkdir -p tmp_entire
seqkit fx2tab /path/to/SQANTI3_output/GRCh38_gencode38_classification.filtered_lite.fasta > entire.txt
tac entire.txt | while read line; do
transcript_ID=$(echo $line | awk '{print $1}')
sequence=$(echo $line | awk '{print $2}')
# echo $transcript_ID
# echo $sequence
window=25
seq_len=$(echo ${#sequence})
# echo ${seq_len}
if [ -e "tmp_entire/${transcript_ID}.txt" ]; then
echo "File already exists"
else
echo -e "isoform\tsequence_length\twindow\tdistance_from_cap\tMFE\tcentroid\tMEA" > tmp_entire/${transcript_ID}.txt
for ((i=1; i<=$seq_len-$window; i=i+1)); do
echo -e "$i / $((seq_len-$window+1))"
tmp_seq=${sequence:$i:$window}
# echo $tmp_seq
tmp_seq_len=$(echo ${#tmp_seq})
# echo $tmp_seq_len
# echo $tmp_seq | RNALfold -L $window
echo $tmp_seq | RNAfold --noPS --MEA -p | cut -d " " -f 2- > tmp_ent_res.txt
#first line is just the seq - skip
#second line has the MFE
MFE=$( sed -n "2p" tmp_ent_res.txt | sed -e "s/(//g" | sed -e "s/^ //g" | awk '{print $1}' | sed -e "s/)//g" )
#third line has the ensemble energy - discarding for now
#fourth line has the centroid energy
centroid=$( sed -n "4p" tmp_ent_res.txt | sed -e "s/{//g" | sed -e "s/^ //g" | awk '{print $1}' )
#fifth line has the MEA
MEA=$( sed -n "5p" tmp_ent_res.txt | sed -e "s/{//g" | sed -e "s/^ //g" | awk '{print $1}' )
echo -e "${transcript_ID}\t${seq_len}\t${window}\t$((i-1))\t${MFE}\t${centroid}\t${MEA}" >> tmp_entire/${transcript_ID}.txt
done || true
rm tmp_ent_res.txt
echo -e 'finish!'
fi
done

cp tmp_entire/* /path/to/RNA_folding_energy/
