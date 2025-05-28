#!/bin/sh
#$ -S /bin/sh

# IUPred3
## https://iupred2a.elte.hu/help_new

tool="SQANTI3"
path_to_aafasta="/path/to/GRCh38_gencode38_classification.filtered_lite_aa.fasta"

grep "^>" ${path_to_aafasta} | sed -e "s/^>//g" | awk '{print $1}' > read_id_${tool}.list
mkdir --parents /path/to/IUPred3_output/${tool}
# extract aa sequence of interested read
cat read_id_${tool}.list | while read read_id; do

echo -e "${tool}, ${read_id}"
seqkit grep -nrp ${read_id} ${path_to_aafasta} > /path/to/tmp/${read_id}_${tool}.fa

/path/to/bio/python3/bin/python3 /path/to/tools/iupred3/iupred3.py \
--anchor \
/path/to/tmp/${read_id}_${tool}.fa \
long \
> /path/to/IUPred3_output/${tool}/${read_id}.txt
done

if [ -e "/path/to/IUPred3_output/${tool}_summary.txt" ]; then
rm /path/to/IUPred3_output/${tool}_summary.txt
fi
touch /path/to/IUPred3_output/${tool}_summary.txt
echo -e "isoform\tidr_length\tanchor2_length" >> /path/to/IUPred3_output/${tool}_summary.txt
cat read_id_${tool}.list | while read read_id; do
total=$(wc -l read_id_${tool}.list | awk '{print $1}')
if [ -e "/path/to/IUPred3_output/${tool}/${read_id}.txt" ]; then
n=$(grep -w -n ${read_id} read_id_${tool}.list | sed -e "s/:/\t/g" | cut -f 1)
N=$(grep -v "^#" /path/to/IUPred3_output/${tool}/${read_id}.txt | awk '$3 > 0.5 {print $0}' | wc -l)
M=$(grep -v "^#" /path/to/IUPred3_output/${tool}/${read_id}.txt | awk '$4 > 0.5 {print $0}' | wc -l)
echo -e "${n}/${total}, ${read_id}\t${N}\t${M}"
echo -e "${read_id}\t${N}\t${M}" >> /path/to/IUPred3_output/${tool}_summary.txt
fi
done
